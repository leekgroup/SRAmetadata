/* asymptote.cpp
   Abhi Nellore
   September 2, 2015

   Consider an exon-exon junction filter where a junction is kept if and only
   if it is found in >= some proportion K of RNA-seq samples analyzed.
   This program takes S random samples of N RNA-seq samples of input junction
   counts from a total of Q RNA-seq samples at various values of K and finds
   the number of junctions J that make it past the filter, writing a table
   to stdout whose rows have the following format.

   N <tab> K <tab> J

   The expected input from stdin is the junction output of Rail for all of SRA.
   Note that some constants dependent on this input data are baked into the
   code below.

   Requires Boost. Random seed is first and only command-line parameter.
   We used seeds on [0, 49]
   */

#include<iostream>
#include<vector>
#include<bitset>
#include<string>
#include<random>
#include "boost/multi_array.hpp"

const unsigned int INTRON_COUNT = 42882032; // Bake intron count in so bitsets can be preallocated
const unsigned int SAMPLE_COUNT = 21506; // Total number of samples
const unsigned int RANDOM_COUNT = 1; // Number of random samples S to take at each value of N
const unsigned int SAMPLE_INTERVAL = 500; // Minimum value of N defined above as well as interval between successive values pf N studied
const unsigned int SAMPLE_MAX = 21500; // Maximum value of N defined above as well as interval between successive values of N studied
const double PROPORTION_INTERVAL = 0.025; // Minimum value of K defined above as well as interval between successive values of K studied
const double PROPORTION_MAX = 0.075; // Maximum value of K studied

/* Sampling from range without replacement implementation inspired
by http://stackoverflow.com/questions/28287138/c-randomly-sample-k-numbers-from-range-0n-1-n-k-without-replacement */

std::bitset<SAMPLE_COUNT> randomSet(int r, std::default_random_engine& gen)
{
    std::bitset<SAMPLE_COUNT> randomSample;
    for (int k = SAMPLE_COUNT - r; k < SAMPLE_COUNT; ++k) {
        int v = std::uniform_int_distribution<>(1, k)(gen);

        // there are two cases.
        // v is not in candidates ==> add it
        // v is in candidates ==> well, k is definitely not, because
        // this is the first iteration in the loop that we could've
        // picked something that big.
        if (!randomSample[v]) {
            randomSample.set(v);
        } else {
            randomSample.set(k);
        }
    }
    return randomSample;
}

std::bitset<SAMPLE_COUNT> intronFromLine(const std::string &str) {
   std::bitset<SAMPLE_COUNT> intronInSampleQ;
   if (!str.length()) return intronInSampleQ;
   int tabCount = 0, i = 0;
   while (tabCount < 6) {
      if (str[i] == '\t') tabCount++;
      i++;
   }
   int startIndex = i, numLength = 0;
   while (str[i] != '\t') {
      if (str[i] == ',') {
         intronInSampleQ.set(std::stoi(str.substr(startIndex, numLength)));
         startIndex = i + 1;
         numLength = 0;
      } else {
         numLength += 1;
      }
      i++;
   }
   intronInSampleQ.set(std::stoi(str.substr(startIndex, numLength)));
   return intronInSampleQ;
}

int main(int argc, char **argv) {
   const unsigned int SEED = atoi(argv[1]);
   int rowCount = SAMPLE_MAX / SAMPLE_INTERVAL;
   int columnCount = PROPORTION_MAX / PROPORTION_INTERVAL + 1;
   std::cerr << "Reservoir sampling to obtain random bitsets..." << std::endl;
   std::default_random_engine gen(SEED);
   std::vector<std::bitset<SAMPLE_COUNT> > randomSamples(rowCount * RANDOM_COUNT);
   std::vector<int > randomSampleCounts(rowCount * RANDOM_COUNT);
   int currentRandomSampleCount;
   int k = 0;
   for (int i = 0; i < rowCount; i++) {
      currentRandomSampleCount = (i + 1) * SAMPLE_INTERVAL;
      for (int j = 0; j < RANDOM_COUNT; j++) {
         randomSamples[k] = randomSet(currentRandomSampleCount, gen);
         randomSampleCounts[k] = randomSamples[k].count();
         k += 1;
      }
   }
   std::cerr << "Reading junctions and applying filters..." << std::endl;
   typedef boost::multi_array<int, 2> array_type;
   typedef array_type::index index;
   array_type junctionCounts(boost::extents[randomSamples.size()][columnCount]);
   std::fill( junctionCounts.origin(), junctionCounts.origin() + junctionCounts.size(), 0 );
   double proportion;
   int numSamplesWithJunction;
   std::string str;
   while (getline(std::cin, str)) {
      std::bitset<SAMPLE_COUNT> intronInSampleQ = intronFromLine(str);
      for (int k = 0; k < columnCount; k++) {
         proportion = k * PROPORTION_INTERVAL;
         for (int i = 0; i < randomSamples.size(); i++) {
            numSamplesWithJunction = (randomSamples[i] & intronInSampleQ).count();
            if (numSamplesWithJunction >= proportion * randomSampleCounts[i] && numSamplesWithJunction > 0) {
               junctionCounts[i][k]++;
            }
         }
      }
   }
   std::cerr << "Dumping output..." << std::endl;
   for (int k = 0; k < columnCount; k++) {
      proportion = k * PROPORTION_INTERVAL;
      for (int i = 0; i < randomSampleCounts.size(); i++) {
         std::cout << std::to_string(randomSampleCounts[i]) << '\t' << std::to_string(proportion) << '\t'
            << std::to_string(junctionCounts[i][k]) << std::endl;
      }
   }
   std::cerr << "Done." << std::endl;
   return 0;
}
/* asymptote.cpp
   Abhi Nellore
   September 2, 2015

   Consider an exon-exon junction filter where a junction is kept if and only
   if it is found in >= some proportion K of RNA-seq samples analyzed.
   This script takes S random samples of N RNA-seq samples of input junction
   counts from a total of Q RNA-seq samples at various values of K and finds
   the number of junctions J that make it past the filter, writing a table
   to stdout whose rows have the following format.

   N <tab> K <tab> J

   The expected input from stdin is the junction output of Rail for all of SRA.
   Note that some constants dependent on this input data are baked into the
   code below.
   */

#include<iostream>
#include<vector>
#include<bitset>
#include<string>
#include<random>

const unsigned int INTRON_COUNT = 42882032; // Bake intron count in so bitsets can be preallocated
const unsigned int SAMPLE_COUNT = 21506; // Total number of samples
const unsigned int RANDOM_COUNT = 50; // Number of random samples S to take at each value of N
const unsigned int SAMPLE_INTERVAL = 100; // Minimum value of N defined above as well as interval between successive values pf N studied
const unsigned int SAMPLE_MAX = 20000; // Maximum value of N defined above as well as interval between successive values of N studied
const unsigned double PROPORTION_INTERVAL = 0.005; // Minimum value of K defined above as well as interval between successive values of K studied
const unsigned double PROPORTION_MAX = 0.03; // Maximum value of K studied
const unsigned int SEED = 5; // Or whatever; we used 5 for reproducibility

void addIntron(const std::string &str, int columnIndex, std::vector<std::bitset<INTRON_COUNT> > &introns) {
   if (!str.length()) return;
   size_t tabCount = 0, i = 0;
   while (tabCount < 3) {
      if (str[i] == '\t') tabCount++;
      i++;
   }
   int startIndex = i, numLength = 0;
   while (str[i] != '\t') {
      if (str[i] == ',') {
         introns[std::stoi(str.substr(startIndex, numLength))].set(columnIndex);
         startIndex = i + 1;
         numLength = 0;
      } else {
         numLength += 1;
      }
      i++;
   }
   introns[std::stoi(str.substr(startIndex, numLength))].set(columnIndex);
   return;
}

int main() {
   std::cerr << "Initializing data structures..." << std::endl;
   std::vector<std::bitset<INTRON_COUNT> > introns(SAMPLE_COUNT);
   int rowCount = SAMPLE_MAX / SAMPLE_INTERVAL
   int columnCount = PROPORTION_MAX / PROPORTION_INTERVAL
   std::vector<std::double<columnCount> > countMeans(rowCount);
   std::vector<std::double<columnCount> > countStandardDeviations(rowCount);
   std::cerr << "Done." << std::endl;
   std::string str;
   int count = 0;
   std::cerr << "Loading junctions..." << std::endl;
   while (getline(std::cin, str)) {
      addIntron(str, count, introns);
      count++;
   }
   std::cerr << "Done." << std::endl;
   std::cerr << "Counting junctions..." << std::endl;
   std::default_random_engine generator(SEED);
   std::uniform_int_distribution<int> distribution(0, SAMPLE_COUNT - 1);
   return 0;
}
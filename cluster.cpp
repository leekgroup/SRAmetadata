/* Performs correlation clustering on "itn" output of Rail
   Abhi Nellore
   June 17, 2015 */
#include<iostream>
#include<vector>
#include<bitset>
#include<string>
#include<random>

const int INTRON_COUNT = 11898514; // Bake intron count in so bitsets can be preallocated
const float THRESHOLD = 0.8; // Jaccard index >= this value is + edge; else - edge
const int SAMPLE_COUNT = 3000;
const unsigned int SEED = 5;

void addIntron(const std::string &str, int columnIndex, std::vector<std::bitset<INTRON_COUNT> > &introns) {
   if (!str.length()) return;
   int tabCount = 0, i = 0;
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
   return;
}

int main() {
   std::cerr << "Initializing data structures..." << std::endl;
   std::vector<std::bitset<INTRON_COUNT> > introns(SAMPLE_COUNT);
   std::cerr << "Done." << std::endl;
   std::string str;
   int count = 0;
   std::cerr << "Loading introns..." << std::endl;
   while (getline(std::cin, str)) {
      addIntron(str, count, introns);
      count++;
   }
   std::cerr << "Done." << std::endl;
   std::cerr << "Clustering..." << std::endl;
   std::default_random_engine generator(SEED);
   std::uniform_int_distribution<int> distribution(0, SAMPLE_COUNT - 1);
   std::vector<int> unclustered(SAMPLE_COUNT);
   std::iota (std::begin(unclustered), std::end(unclustered), 0);
   std::vector<int> newUnclustered(0);
   std::vector<int> newCluster(0);
   int pivot, unioned, intersected;
   float jaccard;
   while (unclustered.size() > 0) {
      pivot = distribution(generator);
      if (!introns[pivot].count()) {
         // No introns? x it out
         std::cout << "x " << std::to_string(pivot) << std::endl;
         continue;
      }
      newCluster.push_back(pivot);
      for (auto &i : unclustered) {
         if (pivot == i) continue;
         intersected = (introns[pivot] & introns[i]).count();
         unioned = (introns[pivot] | introns[i]).count();
         if (unioned) {
            jaccard = (float) intersected / (float) unioned;
         } else {
            jaccard = 0.0;
         }
         if (jaccard < THRESHOLD) {
            newUnclustered.push_back(i);
         } else {
            newCluster.push_back(i);
         }
      }
      unclustered.swap(newUnclustered);
      newUnclustered.clear();
      for (auto &i : newCluster) {
         std::cout << std::to_string(i) << ' ';
      }
      std::cout << std::endl;
      newCluster.clear();
   }
   return 0;
}
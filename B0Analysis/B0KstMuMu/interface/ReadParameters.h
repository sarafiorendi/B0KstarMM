#ifndef _READPARAMETERS
#define _READPARAMETERS

#include <fstream>
#include <string>
#include <vector>


class ReadParameters
{
  
  
 public:
  ReadParameters (const char* InputFileName, std::string howOpen = "in");
  ~ReadParameters () {fileHandler.close();};
  void ReadFromFile (unsigned int blockNumber, std::vector<std::string>* ParVector);
  void SaveToFile   (unsigned int blockNumber, std::vector<std::string>* ParVector);


 private:
  std::fstream fileHandler;
  
};

#endif

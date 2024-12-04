#include<random>
#include"FitMagneticArea.h"

inline std::mt19937& getRngb() {
  static std::random_device rd;
  static std::mt19937 rng(rd());
  return rng;
}

inline double generateRandomDoubleb(double first, double last) {
  std::uniform_real_distribution<double> dist(first, last);
  return dist(getRngb());
}

void readAreas(const std::string& filename, std::vector<FieldParameters>& fb, std::vector<double>& areas) {

  double mu0       = 4*M_PI*1e-7*1e9;
  std::ifstream infile(filename);
  if (!infile) {
    throw std::runtime_error("Could not open file " + filename);
  }
  
  std::string line;
  while (std::getline(infile, line)) {
    // Ignore comment lines starting with #
    if (line.empty() || line[0] == '#') {
      continue;
    }
    
    std::istringstream iss(line);
    double freq, field, area;
    if (!(iss >> freq >> field >> area)) {
      throw std::runtime_error("Error reading line: " + line);
    }
    
    FieldParameters fp; fp.freq = freq*1e-6; fp.field = field*1e-6*mu0;
    fb.push_back(fp);
    areas.push_back(area/(field));
  }
  
  infile.close();
}

std::vector<std::map<std::string, double>> generateInitialGuesses(int nTemperatures){

  //Be careful, here the standard deviation is in units of rc: std = rc*std_in
  std::vector<std::map<std::string, double>> initialGuesses(nTemperatures);
  auto now    = std::chrono::system_clock::now();
  auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
  
  fori(0, nTemperatures){
    double rc            = generateRandomDoubleb(5,50);
    double msat          = generateRandomDoubleb(0.00005,0.0005);
    double rcstd         = generateRandomDoubleb(0.05,1);
    double numParticles  = generateRandomDoubleb(0.01,10);
    std::map<std::string, double> initialGuess = {
      {"coreRadius", rc},
      {"msat", msat},
      {"numParticles", numParticles},
    };
    initialGuesses[i] = initialGuess;    
  }
  return initialGuesses;
}




int main(int argc, char* argv[]){
  
  std::vector<double> MCtemperatures = {1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1};
  std::vector<double> jumpSize(MCtemperatures.size(), 0.1);
  auto initialGuesses = generateInitialGuesses(MCtemperatures.size());

  ParallelTempering::Parameters ptPar;
  ptPar.temperatures   = MCtemperatures;
  ptPar.numStepsFinish = 2000;
  ptPar.numStepsSwap   = 100;
  ptPar.maxIterations  = 2000;
  ptPar.jumpSize       = jumpSize;
  ptPar.tolerance      = 0.00001;
  ptPar.numThreads     = 5;
  ptPar.printSteps     = 1000;
  
  
  GaussNewton::Parameters gnPar;
  gnPar.maxIterations     = 500;
  gnPar.tolerance         = 1e-4;
  gnPar.regularization    = 1e-6;
  gnPar.printSteps        = 5;
  
  std::string filename = std::string(argv[1]);
  std::vector<FieldParameters> fb;
  std::vector<double> areas;
  readAreas(filename, fb, areas);

  std::map<std::string, double> extraParameters = {
    {"kBT", 0.0041124},
    {"viscosity", 0.0009},
    {"coatingWidth", 0.0},
  };
  
  auto params = fitAreas(fb, areas, ptPar, gnPar,
                         initialGuesses, extraParameters);
  
  for (const auto& pair : params.parameters) {
    std::cout << pair.first << ": " << pair.second << std::endl;
  }
}

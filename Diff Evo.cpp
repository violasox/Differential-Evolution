// compiler_options_EHA.cpp
// compile with: /EHa

// Diff Evo.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <map>
#include <ctime>
#include <math.h>
#include <random>
#include <iomanip>
#include "rr/rrRoadRunner.h"
#include "rr/rrExecutableModel.h"

using std::vector;
using std::cout;
using std::endl;
using std::string;

class Member {
public:
	std::vector<double> vector;
	double fitness;
	Member(std::vector<double>, double);
	Member();
};

Member::Member(std::vector<double> ve, double fi) {
	vector = ve;
	fitness = fi;
}

Member::Member() {

}


class ParameterInfo {
public:
	int index;
	bool fixed;
	double minimum;
	double maximum;
	double log_range[2];
	ParameterInfo(int, bool, double, double, double, double);
	ParameterInfo();
};

ParameterInfo::ParameterInfo(int ind, bool fix = false, double min = 0., double max = 10., double lgmin = 0., double lgmax = 1.) {
	index = ind;
	fixed = fix;
	minimum = min;
	maximum = max;
	log_range[0] = lgmin;
	log_range[1] = lgmax;
}

ParameterInfo::ParameterInfo() {
}

class Optimization {
public:
	rr::ExecutableModel *mdl;
	rr::RoadRunner *roadr;
	std::vector<std::vector<double>> measurements;
	std::vector<std::vector<Member>> islands;
	std::vector<std::string> headings;
	std::vector<std::string> parameters;
	std::vector<double> parameterValues;
	std::map <std::string, ParameterInfo> parameter_map;
	int MaximumGenerations;
	double FitnessThreshold;
	const static int NumberOfIslands = 1;
	int PopulationSize;

	double duration;
	double timeBegin;
	double timeStep;

	Optimization(std::string);
	int LoadMeasurements(string, bool, int);
	int CreateParameterMap();
	int FixParameter(std::string);
	int UnfixParameter(std::string);
	int AssignVectorToModel(std::vector<double>);
	std::vector<Member> CreateIsland();
	int SortIslandsByFitness();
	Member CreateRandomMember();
	Member ChallengeMember(Member, std::vector<Member>, double, double); 
	double GetFitness(std::vector<double>);
	int Run();
};

Optimization::Optimization(std::string modelName) {
	rr::RoadRunner r(modelName);
	roadr = &r;
	rr::ExecutableModel *model = r.getModel();
	mdl = model;
}

// Trying to debug this method, it's a mess
int Optimization::LoadMeasurements(string filename, bool header = true, int decimal_precision = 2) {
	// Read file into data. File has only headings and rows of data with each entry separated by a comma.
	std::ifstream datafile(filename);
	if (datafile.is_open()) { cout << "it's open" << endl; }
	// std::vector<vector<double>> data;
	string line;
	std::vector<std::string> headings;
	std::istream_iterator<std::string> start(datafile), end;
	std::vector<std::string> test(start, end);

	cout << test.size() << endl;

	//bool del = true;
	//bool deldel = true;
	//int count = 0;
	//while (std::getline(datafile, line)) {
	//	if (del) {
	//		del = false;
	//		if (header) {
	//			//std::replace(line.begin(), line.end(), ",", " ");
	//			std::stringstream linestream(line);
	//			std::string value;
	//			while (linestream >> value) {
	//				headings.push_back(value);
	//			}
	//		}
	//		else {
	//			//std::replace(line.begin(), line.end(), ",", " ");
	//			std::stringstream linestream(line);
	//			double value;
	//			while (linestream >> value) {
	//				headings.push_back("Measurement " + std::to_string(value));
	//			}
	//		}
	//		continue;
	//	}
	//	else {
	//		//if (deldel && header) {
	//		//	deldel = false;
	//		//	continue;
	//		//}
	//		cout << "got to the else" << endl;
	//		std::vector<double> miniData;
	//		//std::replace(line.begin(), line.end(), ",", " ");
	//		std::stringstream linestream(line);
	//		double value;
	//		while (linestream >> value) {
	//			miniData.push_back(value);
	//		}
	//		measurements.push_back(miniData);
	//	}
	//}

	cout << "I'm here!" << endl;
	if (measurements.empty()) {
		cout << "no data!" << endl;
	}
	for (auto it1 = measurements.begin(); it1 != measurements.end(); ++it1) {
		for (auto it2 = it1->begin(); it2 != it1->end(); ++it2) {
			cout << *it2 << ", ";
		}
		cout << endl;
	}

	if (test.empty()) { cout << "test is empty!" << endl; }
	return 0;
}

int Optimization::CreateParameterMap() {
	std::map <std::string, ParameterInfo> parameter_map;
	int paramSize = mdl->getNumGlobalParameters();
	double *values;
	int *indx;
	values = new double[paramSize];
	indx = new int[paramSize];
	mdl->getGlobalParameterValues(paramSize, indx, values);
	for (int i = 0; i != paramSize; ++i) {
		parameters.push_back(mdl->getGlobalParameterId(i));
		parameterValues.push_back(values[i]);
	}
//	string* parameterString = mdl->getGlobalParameterIds();
	for (int i = 0; i != parameters.size(); ++i) {
		ParameterInfo info = ParameterInfo(parameter_map.size());
		std::string key = parameters[i];
		parameter_map[key] = info;
	}
	return 0;
}

int Optimization::FixParameter(std::string key) {
	parameter_map[key].fixed = true;
	return 0;
}

int Optimization::UnfixParameter(std::string key) {
	parameter_map[key].fixed = false;
	return 0;
}

int Optimization::AssignVectorToModel(std::vector<double> ve) {
	return 0;
}

std::vector<Member> Optimization::CreateIsland() {
	std::vector<Member> island;
	for (int i = 0; i != PopulationSize; ++i) {
		island.push_back(CreateRandomMember());
	}
	return island;
}

int Optimization::SortIslandsByFitness() {
	std::vector<std::map<double, Member>> islandFitness;
	for (int i = 0; i != islands.size(); ++i) {
		std::map<double, Member> islandMembers;
		for (int j = 0; j != islands[i].size(); ++j) {
			islandMembers[islands[i][j].fitness] = islands[i][j];
		}
		islandFitness.push_back(islandMembers);
	}
	
	islands.clear();
	for (int i = 0; i != islandFitness.size(); ++i) {
		std::vector<Member> island;
		for (auto iter = islandFitness[i].begin(); iter != islandFitness[i].end(); ++iter) {
			island.push_back(iter->second);
		}
		islands.push_back(island);
	}
	return 0;
}

Member Optimization::CreateRandomMember() {
	bool success = false;
	std::vector<double> vec;
	double fitness = 0;
	while (!success) {
		try {
			std::vector<double> paramValues;
			int paramSize = mdl->getNumGlobalParameters();
			double *values;
			int *indx;
			values = new double[paramSize];
			indx = new int[paramSize];
			mdl->getGlobalParameterValues(paramSize, indx, values);
			for (int i = 0; i != paramSize; ++i) {
				parameterValues.push_back(values[i]);
			}

			// string* parameterString = mdl->getGlobalParameterIds();
			// int index = 0;
			std::random_device rd; // obtain a random number from hardware
			std::mt19937 eng(rd()); // seed the generator

			for (int i = 0; i != parameters.size(); ++i) {
				if (parameter_map[parameters[i]].fixed) {
					// int ind = std::distance(mdl->getGlobalParameterIds(), std::find(mdl->getGlobalParameterIds(), *parameterString));
					paramValues[i] = parameterValues[i];
				}
				else {
					int rangeMin = parameter_map[parameters[i]].log_range[0];
					int rangeMax = parameter_map[parameters[i]].log_range[1];
					std::uniform_int_distribution<> distr(rangeMin, rangeMax); // define the range
					paramValues[i] = std::pow(10, distr(eng));
				}
			}
			fitness = GetFitness(paramValues);
			if (paramValues.size() != parameters.size()) {
				throw "Error creating the initial optimization vector. Not all parameters were assigned a value.";
			}
			success = true;
			cout << "Created a random member.\tFitness: " << fitness << "\tVector: " << endl;
			std::vector<double>::iterator iter = paramValues.begin();
			while (iter != paramValues.end()) {
				cout << *iter << endl;
				iter++;
			}
		}
		catch (std::exception err) {
			cout << err.what() << endl;
		}
	}

	return Member(vec, fitness);
}

Member Optimization::ChallengeMember(Member original, std::vector<Member> samples, double CR = 0.6, double F = 0.8) {
	std::vector<double> o = original.vector;
	std::vector<double> a = samples[0].vector;
	std::vector<double> b = samples[1].vector;
	std::vector<double> c = samples[2].vector;
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 eng(rd()); // seed the generator

	// Basic mutation
	std::vector<double> newVector;
	for (int i = 0; i != o.size(); ++i) {
		double randDouble = std::rand() / RAND_MAX;
		if (randDouble <= CR) {
			double v = a[i] + F * (b[i] + c[i]);
			if (v > 0) {
				newVector.push_back(v);
			}
			else {
				newVector.push_back(o[i] / 2.);
			}
		}
		else {
			newVector.push_back(o[i]);
		}
		// Check that all vector values conform to the parameter constraints.
		// std::string* parameterString = mdl->getGlobalParameterIds();
		// Reset to original if value was supposed to be fixed.
		if (parameter_map[parameters[i]].fixed) {
			newVector[i] = o[i];
		}
		// Resample within log_range if value exceeds bounds.
		if (newVector[i] < parameter_map[parameters[i]].minimum || newVector[i] > parameter_map[parameters[i]].maximum) {
			int rangeMin = parameter_map[parameters[i]].log_range[0];
			int rangeMax = parameter_map[parameters[i]].log_range[1];
			std::uniform_int_distribution<> distr(rangeMin, rangeMax); // define the range
			newVector[i] = std::pow(10, distr(eng));
		}

		try {
			double fitness = GetFitness(newVector);
			if (fitness < original.fitness) {
				return Member(newVector, fitness);
			}
			else { return original; }
		}
		catch (std::exception) { 
			return original; 
		}
	}
	return original;
}

double Optimization::GetFitness(std::vector<double> ve) {
	mdl->reset();
	AssignVectorToModel(ve);

	rr::BasicDictionary options;
	options.setItem("start", timeBegin);
	options.setItem("duration", duration);
	options.setItem("steps", timeStep);
	//int numRows = roadr->simulate(&options)->numRows;
	//int numCols = roadr->simulate(&options)->numCols;
	const ls::DoubleMatrix *simulatedResults = roadr->simulate();
	int cols = simulatedResults->numCols();
	int rows = simulatedResults->numRows();

	double fitness = 0;
	for (int i = 0; i != rows; ++i) {
		for (int j = 1; j != cols; ++j) {
			double observed = (*simulatedResults)(i, j);
			double expected = measurements[i][j];
			fitness += std::pow((observed - expected), 2);
		}
	}

	return fitness;
}

int Optimization::Run() {
	cout << "\n\nStarting differential evolution routine..." << endl;
	time_t start_time = std::time(0);
	PopulationSize = __max(PopulationSize, 12);
	for (int i = 0; i != NumberOfIslands; ++i) {
		islands.push_back(CreateIsland());
	}

	int generation_count = 0;
	while (generation_count < MaximumGenerations) {
		++generation_count;
		std::vector<std::vector<Member>> populationSamples;
		for (int i = 0; i != islands.size(); ++i) {
			for (int j = 0; j != islands[i].size(); ++j) {
				std::vector<Member> sample;
				double numLeft = 3;
				double size = islands[i].size();

				for (int h = 0; h != islands[i].size(); ++h) {
					double probability = numLeft / size;
					double result = std::rand() / RAND_MAX;
					if (result < probability) {
						sample.push_back(islands[i][h]);
						--numLeft;
					}
					--size;
				}
				populationSamples.push_back(sample);
			}
			for (int index = 0; index != islands[i].size(); ++index) {
				islands[i][index] = ChallengeMember(islands[i][index], populationSamples[index]);
			}
		}
		SortIslandsByFitness();
		if (generation_count % 10 == 0) {
			cout << endl;
			for (int i = 0; i != islands.size(); ++i) {
				for (int j = 0; j != islands[i].size(); ++j) {
					cout << "Fitness: " << std::setprecision(4) << islands[i][j].fitness << "\tVector: ";
					for (auto iter = islands[i][j].vector.begin(); iter != islands[i][j].vector.end(); ++iter) { cout << *iter; }
					cout << endl;
				}
			}
		}
		std::vector<double> fitnessVector;
		for (int i = 0; i != islands.size(); ++i) {
			fitnessVector.push_back(islands[i][0].fitness);
		}
		double minimum = *std::min_element(fitnessVector.begin(), fitnessVector.end());
		if (minimum < FitnessThreshold) { break; }
	}
	// Print final report
	for (int i = 0; i != islands.size(); ++i) {
		for (int j = 0; j != islands[i].size(); ++j) {
			cout << "Fitness: " << std::setprecision(4) << islands[i][j].fitness << "\tVector: ";
			for (auto iter = islands[i][j].vector.begin(); iter != islands[i][j].vector.end(); ++iter) { cout << *iter; }
			cout << endl;
		}
	}
	cout << "\nCompleted differential evolution routine. Total time: " << std::time(0) - start_time;
	return 0;
}



using namespace rr;
int main()
{
	try {
		Optimization opt("C:/Users/user/Documents/tellurium-winpython/Model Fitting/Model 4.xml");

		opt.LoadMeasurements("C:/Users/user/Documents/tellurium-winpython/Model Fitting/experimental_data_test3.dat");
	}

	catch (std::exception ex) {
		cout << ex.what();
	}

	//	opt.MaximumGenerations = 200;
	//	opt.FitnessThreshold = 1e-3;
	//	opt.NumberOfIslands = 1;
	//	opt.PopulationSize = 12;
	//	opt.MigrationFrequency = 0.3;
	//	opt.NumberOfMigrants = 1;
	//	opt.SelectionGroupSize = 3;
	//	opt.ReplacementGroupSize = 3;
	//	opt.MigrationTopology[0] = 0;

	//opt.LoadMeasurements("file.txt");
	//opt.CreateParameterMap();

	/*
		// Create exact results (still in Python)
		result = opt.model.simulate(0, 15, 301)
		result = np.asarray([[y for y in x] for x in result])
		entities = [x.replace('[', '').replace(']', '') for x in opt.model.selections]
		opt.model.plot()
		print('{}\t{}\t{}\t{}\t'.format(*entities))
		for i in range(len(result)) :
			print('{}\t{}\t{}\t{}\t'.format(*result[i]))
	
	*/

	//opt.Run();

}

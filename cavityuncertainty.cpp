#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <random>
#include <iostream>

#include "src/Structs.hpp"

/*
 * http://stackoverflow.com/questions/236129/split-a-string-in-c
 */
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

/*
 * http://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
 */
template <typename T, typename Compare>
std::vector<int> sort_permutation(
	std::vector<T> const& vec,
	Compare compare)
{
	std::vector<int> p(vec.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](int i, int j){ return compare(vec[i], vec[j]); });
	return p;
}

template <typename T>
std::vector<T> apply_permutation(
	std::vector<T> const& vec,
	std::vector<int> const& p)
{
	std::vector<T> sorted_vec(p.size());
	std::transform(p.begin(), p.end(), sorted_vec.begin(),
		[&](int i){ return vec[i]; });
	return sorted_vec;
}

void Expectation(std::vector<Real>& expec, std::vector<Real> values)
{
	Real sum = 0.0;
	Real n = values.size();
	for (auto x : values)
	{
		sum += x;
	}
	sum /= n;
	expec.push_back(sum);
}

void Variance(std::vector<Real>& variance, Real expec, std::vector<Real> values)
{
	Real sum = 0.0;
	for (auto x : values)
	{
		sum += pow(x - expec, 2.0);
	}
	variance.push_back(sum);
}

void CreateGnuplotOutput(std::vector<Real>* expec, std::vector<Real>* variance, std::vector<Real>& times)
{
	// time E1 E2 E3 E4 E5 E6 V1 V2 V3 V4 V5 V6
	std::string line = "";
	for (unsigned int i = 0; i < times.size(); i++)
	{
		line.append(std::to_string(times[i]));
		line.append(" ");
		for (unsigned int j = 0; j < 6; j++)
		{
			line.append(std::to_string(expec[j][i]));
			line.append(" ");
		}
		for (unsigned int j = 0; j < 6; j++)
		{
			line.append(std::to_string(variance[j][i]));
			if(j < 5) line.append(" ");
			else line.append("\n");
		}
	}
	std::ofstream output("plot.data");
	if (output.is_open())
	{
		output << line;
		output.close();
	}
}

#if defined(__linux)
#include "sys/stat.h"
#endif
#if defined(_WIN64)
#include <Windows.h>
#endif
#if defined(_WIN32)
#include <Windows.h>
#endif
void checkOutputPath()
{
	std::string filename("./uncertainty/");

	// Test if the directory exits, if not create it.
	std::filebuf test_dir;
	test_dir.open(const_cast < char *>(filename.c_str()), std::ios::out);
	if (!test_dir.is_open())
	{
		// Directory doesn't exist.
#if defined(_WIN64)
		CreateDirectory(filename.c_str(), NULL);
#elif defined(_WIN32)
		CreateDirectory(filename.c_str(), NULL);
#elif defined(__linux)
		mkdir(filename.c_str(), 0700);
#endif
	}
}


int main()
{
	checkOutputPath();

	/**
	 * Use max_sim to control the number of simulations that are started.
	 * First we create the configuration files for the simulations, by
	 * choosing a random re number. This is handled differently for linux
	 * and windows
	 */
	Real min_re = 1000;
	Real max_re = 2000;
	Real re = 1000;
	unsigned int max_sim = 1500;
	std::string params;

/* ATTENTION: FROM HERE UNTIL THE OTHER 'ATTENTION' POINT,
 * YOU MAY WANT TO RUN THE CODE PARTS (THIS AND THE CODE AFTER THIS)
 * SEPARATELY. 
 * SO: FIRST GENERATE THE CONFIG FILES, RUN SIMULATIONS ON A LOT OF COMPUTERS,
 *     THEN LET THE SECOND PART OF THIS PROGRAM COMPUTE EXPECTED VALUES AND VARIANCE
 */
	std::string re_numbers("");

	/**
	 * First delete the old Data from previous simulations.
	 */
	std::string filename_rm;
	for (unsigned int i = 0; i < max_sim; i++)
	{

		filename_rm = "./uncertainty/uncertainty";
		filename_rm.append(std::to_string(i));
		filename_rm.append(".txt");
		remove(filename_rm.c_str());
	}

	std::default_random_engine generator;
	std::random_device rd;
	std::normal_distribution<Real> distribution(1500, (max_re - min_re)/6.0);

	for (unsigned int i = 0; i < max_sim; i++)
	{
		re = distribution(rd);
		re_numbers.append(std::to_string(re) + " \n");
		
#if defined(__linux)
		params = "./bin/cavitybaker 0 --re=";
		params.append(std::to_string(re));
		params.append(" --iterMax=100 --tEnd=16.5 --xLength=1 --yLength=1 --iMax=128 --jMax=128 --name=uncertainty");
		params.append(std::to_string(i));
#endif
#if defined(_WIN64)
		params = "cavitybaker.exe 0 --re=";
		params.append(std::to_string(re));
		params.append(" --iterMax=100 --tEnd=16.5 --xLength=1 --yLength=1 --iMax=128 --jMax=128 --name=uncertainty");
		params.append(std::to_string(i));
#endif
#if defined(_WIN32)
		params = "cavitybaker.exe 0 --re=";
		params.append(std::to_string(re));
		params.append(" --iterMax=100 --tEnd=16.5 --xLength=1 --yLength=1 --iMax=128 --jMax=128 --name=uncertainty");
		params.append(std::to_string(i));
#endif

		std::system(params.c_str());
	}

	std::cout << re_numbers;
/* ATTENTION: FIRST CODE PART ENDS HERE. */
	//return 0;

/* ATTENTION: YOU MAY WANT TO DO THE FOLLOWING COMPUTATIONS
 * IN A DISTRIBUTED MANNER. */
	/**
	 * After the configuration files have been created the simualtions are
	 * done sequentially. Each simulation outputs the values for u and v
	 * at the positions (120,5),(64,64),(5,120). 
	 */
	for (unsigned int i = 0; i < max_sim; i++)
	{
#if defined(__linux)
		params = "./bin/cavitydriver uncertainty";
		params.append(std::to_string(i));
#endif
#if defined(_WIN64)
		params = "cavitydriver.exe uncertainty";
		params.append(std::to_string(i));
#endif
#if defined(_WIN32)
		params = "cavitydriver.exe uncertainty";
		params.append(std::to_string(i));
#endif
		std::system(params.c_str());
	}

/* ATTENTION: SECOND PART OF THE CODE STARTS HERE. */
	/**
	 * For each simulation there is one txt file with all values.
	 * In each line there are the values for the three positions at
	 * a certain timestep. Parse this information and store it in
	 * 4 vectors.
	 */
	std::vector<std::pair<Real, Real> > _val_120_5;
	std::vector<std::pair<Real, Real> > _val_64_64;
	std::vector<std::pair<Real, Real> > _val_5_120;
	std::vector<Real> time;
	std::string filename;
	std::string line;
	for (unsigned int i = 0; i < max_sim; i++)
	{
		
		filename = "./uncertainty/uncertainty";
		filename.append(std::to_string(i));
		filename.append(".txt");		
		std::ifstream read_uncer(filename);
		if (read_uncer.is_open())
		{
			while (getline(read_uncer, line))
			{
				std::vector<std::string> elems = split(line, '|');
				std::vector<std::string> val_120_5 = split(elems[0], ',');
				_val_120_5.push_back(std::pair<Real, Real>(atof(val_120_5[0].c_str()), atof(val_120_5[1].c_str())));
				std::vector<std::string> val_64_64 = split(elems[1], ',');
				_val_64_64.push_back(std::pair<Real, Real>(atof(val_64_64[0].c_str()), atof(val_64_64[1].c_str())));
				std::vector<std::string> val_5_120 = split(elems[2], ',');
				_val_5_120.push_back(std::pair<Real, Real>(atof(val_5_120[0].c_str()), atof(val_5_120[1].c_str())));
				time.push_back(atof(elems[3].c_str())); 
			}
			read_uncer.close();
		}
	}

	/**
	 * Now we want to sort this information, that it can be plotted.
	 * First we sort the time vector and store what values have been
	 * swaped. These swaps are than apllied to each of the other three
	 * vectors.
	 */
	auto p = sort_permutation(time,
		[](Real const& a, Real const& b){ return a < b; });

	time = apply_permutation(time, p);
	_val_120_5 = apply_permutation(_val_120_5, p);
	_val_64_64 = apply_permutation(_val_64_64, p);
	_val_5_120 = apply_permutation(_val_5_120, p);

	std::vector<Real> expec[6];
	std::vector<Real> variance[6];
	std::vector<Real> new_time;
	std::vector<Real> values[6];
	Real timestep = time[0];
	new_time.push_back(timestep);
	for (unsigned int i = 0; i < time.size(); i++)
	{
		if (time[i] > timestep)
		{
			for (unsigned int j = 0; j < 6; j++)
			{
				Expectation(expec[j], values[j]);
				Variance(variance[j], expec[j].back(), values[j]);
				values[j].clear();
			}
			timestep = time[i];
			new_time.push_back(timestep);
		}
		values[0].push_back(_val_120_5[i].first);
		values[1].push_back(_val_120_5[i].second);
		values[2].push_back(_val_64_64[i].first);
		values[3].push_back(_val_64_64[i].second);
		values[4].push_back(_val_5_120[i].first);
		values[5].push_back(_val_5_120[i].second);
	}
	for (unsigned int j = 0; j < 6; j++)
	{
		Expectation(expec[j], values[j]);
		Variance(variance[j], expec[j].back(), values[j]);
		values[j].clear();
	}

	CreateGnuplotOutput(expec, variance, new_time);
}

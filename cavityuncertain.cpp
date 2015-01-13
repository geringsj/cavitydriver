#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>

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

int main()
{
	/**
	 * Use max_sim to control the number of simulations that are started.
	 * First we create the configuration files for the simulations, by
	 * choosing a random re number. This is handled differently for linux
	 * and windows (TODO: look at the linux stuff...).
	 */
	int min_re = 1000;
	int max_re = 2000;
	int re = 1000;
	unsigned int max_sim = 2;
	std::string params;
	for (unsigned int i = 0; i < max_sim; i++)
	{
		re = min_re + (rand() % (int)(max_re - min_re + 1));
		
#if defined(__linux)
		params = "./cavitybaker 0 --re=";
		params.append(std::to_string(re));
		params.append(" --iterMax=100 --tEnd=16.5 --xLength=1 --yLength=1 --iMax=128 --jMax=128 --name=uncertainty");
		params.append(std::to_string(i));
#endif
#if defined(_WIN64)
		params = "cavitybaker.exe 0 --re=";
		params.append(std::to_string(re));
		params.append(" --iterMax=100 --tEnd=1.5 --xLength=1 --yLength=1 --iMax=128 --jMax=128 --name=uncertainty");
		params.append(std::to_string(i));
#endif
#if defined(_WIN32)
		params = "cavitybaker.exe 0 --re=";
		params.append(std::to_string(re));
		params.append(" --iterMax=100 --tEnd=1.5 --xLength=1 --yLength=1 --iMax=128 --jMax=128 --name=uncertainty");
		params.append(std::to_string(i));
#endif

		std::system(params.c_str());
	}

	/**
	 * After the configuration files have been created the simualtions are
	 * done sequentially. Each simulation outputs the values for u and v
	 * at the positions (120,5),(64,64),(5,120). 
	 */
	for (unsigned int i = 0; i < max_sim; i++)
	{
#if defined(__linux)
		params = "./cavitydriver uncertainty";
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
}
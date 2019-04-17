#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdint>
#include <unordered_map>
#include "bytell_hash_map.hpp"
#include <random>

typedef std::mt19937 RNG_generator;

size_t hashPair(const std::pair<int, int> & point);

struct customHash /*customHash takes std::pair as input and returns the output of hashPair!!!!*/
{
	public:
		size_t operator()(const std::pair<int,int> & input) const
		{
			return hashPair(input);
		}
}; 

void printPoint(const std::pair<int, int> & p); /* prints point given in pair and data in array*/

void printState(const std::array<bool,4> & arr);

void printLengths(const std::vector<int> & vec);

void printAvailable(const std::vector<std::pair<int,int>> & vec);

void printBurst(const std::vector<std::array<int,3>> & vec);

void printCluster(const ska::bytell_hash_map<std::pair<int, int>, std::array<bool, 4>, customHash> & cluster);

/* void printCluster(const std::unordered_map<std::pair<int, int>, std::array<bool, 4>, customHash> & cluster); */

uint32_t randomChoice(const std::vector<std::pair<int,int>> & vec, RNG_generator & rng); /*returns a random integer in range */

std::pair<int,int> choice(std::vector<std::pair<int,int>> & vec, RNG_generator & rng); /*returns and deletes a random pair from a vector*/

double prob(const std::pair<int, int> & p, const float & p0); /*returns breaking probability*/


std::pair<int, int> getEndSite(const std::pair<int, int> & point, int index); /*gives the end site when pair<x,y> has a bond break on index i=[0,1,2,3] where index s like on a clock. ie 0 is up, 1 right etc*/

int main(int argc, char *argv[])
{
	int t; /* initialize stuff i need.*/
	uint32_t seed_val = time(0); /*seed should not be 2 in general...*/	
	bool full=0;
	const int n = atoi(argv[1]); 
	const float prob_zero = atof(argv[2]);
	const int seed_factor = atoi(argv[3]);
	if (argc == 5)
	{
	std::string fullarg = argv[4]; 
		if (fullarg=="-full")
		{
			full=1;
		}
	}	
	int bursts = 0;

	std::cout << n << " " << prob_zero << " " << full << "\n";

	RNG_generator rng; /* init rng machine*/
	rng.seed(seed_val*seed_factor); 
	std::uniform_real_distribution<double> uniform(0.0, 1.0); /* this we can call to generate random variable on the unit interval */

	std::vector<int> burst_time(n,0); 
	
	/* i think these two can be optimized away... First, find a method that gets me a random key from a map. Second, instead of keeping a burst_list, just dump burst info once a bust is complete. This way i dont have to reallocate memory all the time! This is what makes it FAST^TM*/
	std::vector<std::pair<int,int>> site_list; /* contains all sites */

	std::vector<int> burst_lengths(n,0); /* pre-allocates a vector initally filled w/ zeros of length n */

	ska::bytell_hash_map<std::pair<int,int>, std::array<bool,4>, customHash> cluster_map;
	
	cluster_map[std::pair<int,int> (0,0)] = {0,0,0,0}; /* set (0,0) to empty, and initialize site_list */
	site_list = {std::pair<int,int> (0,0)};

	while (bursts < n)
	{
		/* we choose a random point in site_list to be the start of our burst, and then inits vectors to keep track of the
		 * sites we have available and the ones we have broken. We also add to the timestep and init a burst state bool. */
		uint32_t random_index = randomChoice(site_list, rng);
		std::pair<int,int> initial_site = site_list[random_index];
		std::vector<std::pair<int,int>> available_sites = {initial_site};
		std::vector<std::array<int,3>> broken_bonds;
		t++;
		bool burst_state = 0; /* if burst state == 1, we have emptied the available sites.*/
		while (burst_state == 0)
		{
			std::pair<int,int> site = choice(available_sites, rng); /* removes and returns a random element in availbe_sites */
			std::array<bool,4> site_state = cluster_map[site];
			for (int index=0; index<4; index++)
			{
				if (site_state[index] == 0) /* check that the index is not allready broken... */
				{
					std::pair<int,int>  end_site = getEndSite(site, index);
					if (cluster_map.count(end_site) == 0) /* we check loopless.. */
					{
						double p = uniform(rng);
						double p_tresh = prob(site, prob_zero); /* we check if the bond breaks..*/
						if (p<p_tresh)
						{
						cluster_map[site][index] = 1; /* the iterator open_index.cbegin() .. open_index.cend() gives pointers. *index then returns the value pointed to. (i think..). This is to change the bond at index to be broken on site. We also break on the end site: */
						std::array<bool,4> end_site_state = {0,0,0,0};
						int end_site_index = (index+2)%4; /* this translates [0,1,2,3] -> [2,3,0,1] ie to the opposite side which is the index on the opposite end of the bond*/
						std::array<int,3> broken_bond = {site.first, site.second, index};
						end_site_state[end_site_index] = 1; 
						cluster_map[end_site] = end_site_state; 
						available_sites.push_back(end_site);
						broken_bonds.push_back(broken_bond);
						site_list.push_back(end_site);
						}
					}
				}
			}
			if (available_sites.size() == 0)
			{
				burst_state = 1;
			}
		}
		if (broken_bonds.size() != 0)
		{
			burst_lengths[bursts] = broken_bonds.size();
			burst_time[bursts]=t;
			bursts++;
			if (full==1)
			{
				printBurst(broken_bonds);
			}
		}
	}		
	printLengths(burst_lengths);
	std::cout << std::endl;
}


size_t hashPair(const std::pair<int, int> & point) 
{
	/*function called hash_pair that takes pair (which is from std:: lib) and 
	 * creates a variable point i can reference. pair.first and pair.second 
	 * references the variables....*/
	return (size_t(uint32_t(point.first)) << 32) | size_t(uint32_t(point.second));
	/* we treate the points as unsigned 32bit ints and then assign them to
	 * size_t which just accepts them???. then we shift one over and ??ands??
	 * them together*/
}

std::pair<int, int> getEndSite(const std::pair<int, int> & point, int index)
{
	std::pair<int,int> out;
	if (index == 0)
	{	out.first = point.first;
		out.second = point.second + 1;
	}

	if (index == 1)
	{
		out.first = point.first + 1;
		out.second = point.second;
	}

	if (index == 2)
	{
		out.first = point.first;
		out.second = point.second - 1;
	}

	if (index == 3)
	{
		out.first = point.first -1;
		out.second = point.second;
	}
	return out; 
}

uint32_t randomChoice(const std::vector<std::pair<int,int>> & vec, RNG_generator & rng) /*returns a random integer in range */
{
	std::size_t size = vec.size()-1;
	std::uniform_int_distribution<uint32_t> rng_dist(0, size); /*init a distribution from 0 to size*/
	uint32_t r = rng_dist(rng); 
	return r;
}

std::pair<int,int> choice(std::vector<std::pair<int,int>> & vec, RNG_generator & rng) /*returns and deletes a random pair from a vector*/
{
	uint32_t index = randomChoice(vec, rng); 
	std::pair<int,int> rand_pair = vec[index]; 
	vec.erase(vec.begin() + index);
	return rand_pair;
}


void printLengths(const std::vector<int> & vec)
{
	for (auto it=vec.cbegin(); it<vec.cend(); it++) 
	{
		std::cout <<*it << "\t";
	}
}

void printPoint(const std::pair<int, int> & p)
{
	std::cout << "(" << p.first << ", " << p.second << ")";
}

void printState(const std::array<bool,4> & arr)
{
	std::cout << "[" << arr[0] << ", " << arr[1] << ", " << arr[2] << ", " << arr[3] << "]";
}

void printBurst(const std::vector<std::array<int,3>> & vec)
{
	int size = vec.size();
	for (int i=0; i<size; i++)
	{
		std::cout << "(" << vec[i][0] << " " << vec[i][1] << " " << vec[i][2] << "), ";
	}	
	std::cout << "\n";
}

void printAvailable(const std::vector<std::pair<int,int>> & vec)
{
	int size = vec.size();
	std::cout << "[";
	for (int i=0; i<size; i++)
	{
		std::cout << "(" << vec[i].first << ", " << vec[i].second << "), ";
	}	
	std::cout << "]";
}

void printCluster(const ska::bytell_hash_map<std::pair<int, int>, std::array<bool, 4>, customHash> & cluster)
{	
	std::pair<int, int> key;
	std::array<bool, 4> value;
	for (auto i = cluster.begin(); i != cluster.end(); ++i )/*i iterates over the cluster key-value pairs. i->first will give me the key. i->second will give me the value. (that is completely insane syntax who made this?????*/	
	{ 
		key = i->first;
		value = i->second;
		printPoint(key);
		std::cout << "\t";
		printState(value);
		std::cout << "\n";
	}
}

double prob(const std::pair<int, int> & p, const float & p0)
{
	return p0;	
}

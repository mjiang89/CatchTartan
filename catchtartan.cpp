#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

typedef double real;
typedef long long bits;

struct Cell {
	int sz;
	int *values;
};

struct Tartan {
	vector<int> tlist;
	vector<vector<int> > blist;
	bool *dimlist;
	set<int> *valuesetlist;
};

// t: timestr
// 0: username
// 1: location
// 2: list of phrases
// 3: list of hashtags
// 4: list of urls
// 5: list of rtnames
// 6: list of atnames
const int seed_dim = 2; // seeding with the "phrase" dimension
// 0:time, 1: user, 2: location, 3: phrase, 4: hashtag, 5: URL
const int bounds[10] = {30,2,2,10,2,2,2,2,2,2};
const double _decay = 0.01;

string workdir, datafile, valuefile, outfile;
int num_threads, num_seeds, sz_seed, max_iter;

int num_dims, num_slices;
int *sz_dims; // [num_dims]
int *sz_behaviors; // [num_slices]
vector<Cell *> *behaviorlist; // [num_slices]
vector<string> *tweetlist; // [num_slices]
vector<string> *valuelist; // [num_dims+1] including time

// Toolbox:
// 1) Fast generate a random integer in [0, table_size)
int Rand(unsigned long long &randseed, int table_size) {
	randseed = randseed * 25214903917 + 11;
	return (randseed >> 16) % table_size;
}
// 2) Integer to string
string itos(int x) {
	string ret;
	stringstream ss;
	ss << x;
	ss >> ret;
	return ret;
}
// 3) Split a string [s] with the [delim] into a vector of size [sz]
vector<string> Split(string s, char delim, int sz) {
	vector<string> ret;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) ret.push_back(item);
	if (sz > 0) while (ret.size() < sz) ret.push_back("");
	return ret;
}
// 4) Comparer used for sorting a map-shape vector by the "value"
bool value_comparer(const pair<int, int> &a, const pair<int, int> &b) {
	return a.second > b.second;
}
// 5) Print Tartan
void PrintTartan(struct Tartan &tartan, string title) {
	int i, sz_t, _t, sz_b, sum_b;
	sz_t = tartan.tlist.size();
	cout << title << ": #slices=" << sz_t << endl;
	cout << title << ": slices={";
	for (_t = 0;_t < sz_t;_t++) {
		cout << " " << tartan.tlist[_t];
	}
	cout << " }" << endl;
	cout << title << ": behaviors={";
	sum_b = 0;	
	for (_t = 0;_t < sz_t;_t++) {
		sz_b = tartan.blist[_t].size();
		sum_b += sz_b;
		cout << " " << sz_b;
	}
	cout << " }" << endl;
	cout << title << ": #behaviors=" << sum_b << endl;
	for (i = 0;i < num_dims;i++) {
		cout << title << ": [" << i << "] (" << tartan.dimlist[i] << ") #values="
			<< tartan.valuesetlist[i].size() << endl;
	}	
}
// 6) Get the number of behaviors from a Tartan
int GetNumBehavior(struct Tartan &tartan) {
	int sz_t, _t, sum_b;
	sz_t = tartan.tlist.size();
	sum_b = 0;	
	for (_t = 0;_t < sz_t;_t++) {
		sum_b += tartan.blist[_t].size();
	}
	return sum_b;
}

void LoadData() {
	int i, j, t, value;
	string line;
	vector<string> arr, arr_dim;
	ifstream fr;

	fr.open(workdir+datafile);
	getline(fr, line);
	arr = Split(line, '|', 0);
	num_dims = arr.size()-1;
	num_slices = stoi(arr[0].substr(1));
	sz_dims = new int[num_dims];
	for (i = 0;i < num_dims;i++) {
		sz_dims[i] = stoi(arr[i+1].substr(1));
	}
	sz_behaviors = new int[num_slices];
	for (t = 0;t < num_slices;t++) {
		sz_behaviors[t] = 0;
	}
	behaviorlist = new vector<Cell *>[num_slices];
	tweetlist = new vector<string>[num_slices];
	while (getline(fr, line)) {
		arr = Split(line, '|', num_dims+1);
		t = stoi(arr[0]);
		Cell *cell = new Cell[num_dims];
		for (i = 0;i < num_dims;i++) {
			if (arr[i+1].length() > 0) {
				arr_dim = Split(arr[i+1], ',', 0);
				cell[i].sz = arr_dim.size();
				cell[i].values = new int[cell[i].sz];
				for (j = 0;j < cell[i].sz;j++) {
					value = stoi(arr_dim[j]);
					cell[i].values[j] = value;
				}
			} else {
				cell[i].sz = 0;	
			}
		}
		behaviorlist[t].push_back(cell);
		getline(fr, line);
		sz_behaviors[t]++;
		tweetlist[t].push_back(line);
	}
	fr.close();

	valuelist = new vector<string>[num_dims+1];
	for (i = 0;i < num_dims+1;i++) {
		fr.open(workdir+valuefile+itos(i)+".txt");
		getline(fr, line);
		while (getline(fr, line)) {
			arr = Split(line, '|', 0);
			valuelist[i].push_back(arr[1]);
		}
		fr.close();
	}
}

void GenerateSeed(struct Tartan &tartan, unsigned long long &randseed) {
	int t, d, b, i, j;
	vector<int>::iterator iter;
	tartan.dimlist = new bool[num_dims];
	tartan.valuesetlist = new set<int>[num_dims];
	t = Rand(randseed, num_slices);
	tartan.tlist.push_back(t);
	vector<int> bs;
	i = seed_dim;	
	tartan.dimlist[i] = true;
	for (d = 0;d < sz_seed;d++) {
		b = Rand(randseed, sz_behaviors[t]);
		iter = find(bs.begin(), bs.end(), b);
		if (iter != bs.end()) continue;
		bs.push_back(b);
		{
			if (behaviorlist[t][b][i].sz > 0) {
				for (j = 0;j < behaviorlist[t][b][i].sz;j++) {
					tartan.valuesetlist[i].insert(behaviorlist[t][b][i].values[j]);
				}
			}
		}

	}
	tartan.blist.push_back(bs);
}

void OutputTartan(struct Tartan &tartan, double timecost, int num_iter, int thread_id, int seed_id) {
	int i, j, sz_t, _t, t, b, sum_b, sum_e, value, count;
	map<int, int> *value2count;
	vector<pair<int, int> > value_count;	
	set<int>::iterator iter_set;
	vector<int>::iterator iter_vector;
	map<int, int>::iterator iter_map;
	vector<pair<int, int> >::iterator iter_vector_pair;

	sz_t = tartan.tlist.size();
	sum_b = 0; sum_e = 0;
	value2count = new map<int, int>[num_dims];
	for (_t = 0;_t < sz_t;_t++) {
		t = tartan.tlist[_t];
		sum_b += tartan.blist[_t].size();
		for (iter_vector = tartan.blist[_t].begin();iter_vector != tartan.blist[_t].end();iter_vector++) {
			b = *iter_vector;
			for (i = 0;i < num_dims;i++) {
				if (tartan.dimlist[i]) {
					for (j = 0;j < behaviorlist[t][b][i].sz;j++) {
						value = behaviorlist[t][b][i].values[j];
						iter_set = tartan.valuesetlist[i].find(value);
						if (iter_set != tartan.valuesetlist[i].end()) {
							value2count[i][value]++;
							sum_e++;
						}
					}
				}
			}
		}
	}
	
	string file = workdir+outfile+itos(thread_id)+"_"+itos(seed_id)+".txt";
	cout << "Output " << file << endl;
	ofstream fw(file);
	fw << "thread_id|seed_id|num_behavior|num_entries|num_iter|timecost\n";
	fw << thread_id << "|" << seed_id << "|" << sum_b << "|" << sum_e << "|" << num_iter << "|" << timecost << "\n";
	fw << "sz_time";
	for (i = 0;i < num_dims;i++) {
		fw << "|sz_dim[" << i << "]";
	}
	fw << "\n";
	fw << sz_t;
	for (i = 0;i < num_dims;i++) {
		if (tartan.dimlist[i]) {
			fw << "|" << tartan.valuesetlist[i].size();
		} else {
			fw << "|" << 0;
		}
	}
	fw << "\n";
	for (_t = 0;_t < sz_t;_t++) {
		value = tartan.tlist[_t];
		count = tartan.blist[_t].size();
		fw << '@' << "t" << "|" << count << "|" << value << "|" << valuelist[0][value] << "\n";
	}
	for (i = 0;i < num_dims;i++) {
		if (tartan.dimlist[i]) {
			value_count.clear();
			for (iter_map = value2count[i].begin();iter_map != value2count[i].end();iter_map++) {
				value_count.push_back(*iter_map);
			}
			sort(value_count.begin(), value_count.end(), value_comparer);
			for (iter_vector_pair = value_count.begin();iter_vector_pair != value_count.end();iter_vector_pair++) {
				value = iter_vector_pair->first;
				count = iter_vector_pair->second;
				fw << '@' << i << "|" << count << "|" << value << "|" << valuelist[i+1][value] << "\n";
			}
		}
	}
	for (_t = 0;_t < sz_t;_t++) {
		t = tartan.tlist[_t];
		for (iter_vector = tartan.blist[_t].begin();iter_vector != tartan.blist[_t].end();iter_vector++) {
			b = *iter_vector;
			fw << '#' << tweetlist[t][b] << "\n";
		}
	}
	fw.close();
}

int UpdateBehavior(struct Tartan &tartan) {
	int sz_t, _t, t, b;
	int i, j, value;
	bool isValidBehavior, isValidCell;
	set<int>::iterator iter_set;
	sz_t = tartan.tlist.size();
	for (_t = 0;_t < sz_t;_t++) {
		t = tartan.tlist[_t];
		tartan.blist[_t].clear();
		for (b = 0;b < sz_behaviors[t];b++) {
			isValidBehavior = true;
			for (i = 0;i < num_dims;i++) {
				if (tartan.dimlist[i]) {
					if (behaviorlist[t][b][i].sz == 0) {
						isValidBehavior = false;
						break;
					}
					isValidCell = false;
					for (j = 0;j < behaviorlist[t][b][i].sz;j++) {
						value = behaviorlist[t][b][i].values[j];
						iter_set = tartan.valuesetlist[i].find(value);
						if (iter_set != tartan.valuesetlist[i].end()) {
							isValidCell = true;
							break;
						}
					}
					if (!isValidCell) {
						isValidBehavior = false;
						break;
					}
				}
			}
			if (isValidBehavior) {
				tartan.blist[_t].push_back(b);
			}
		}
	}
	return 0;
}

int UpdateDimension(struct Tartan &tartan, int dim) {
	int sz_t, _t, t, sz_b, _b, b;
	int i, j, value, sz;
	int count, count_pre, count_valid;
	bool isValidBehavior, isValidCell;
	map<int, int> value2count;
	set<int> valueset;
	set<int>::iterator iter_set;
	map<int, int>::iterator iter_map;

	count = 0; count_pre = 0;
	sz_t = tartan.tlist.size();
	for (_t = 0;_t < sz_t;_t++) {
		t = tartan.tlist[_t];
		sz_b = tartan.blist[_t].size();
		for (_b = 0;_b < sz_b;_b++) {
			b = tartan.blist[_t][_b];
			isValidBehavior = true;
			count_valid = 0;
			for (i = 0;i < num_dims;i++) {
				if (tartan.dimlist[i] && i != dim) {
					if (behaviorlist[t][b][i].sz == 0) {
						isValidBehavior = false;
						break;
					}
					isValidCell = false;
					for (j = 0;j < behaviorlist[t][b][i].sz;j++) {
						value = behaviorlist[t][b][i].values[j];
						iter_set = tartan.valuesetlist[i].find(value);
						if (iter_set != tartan.valuesetlist[i].end()) {
							isValidCell = true;
							count_valid++;
						}
					}
					if (!isValidCell) {
						isValidBehavior = false;
						break;
					}
				}
			}
			if (isValidBehavior) {
				count_pre += count_valid;
				sz = behaviorlist[t][b][dim].sz;
				if (sz > 0) {
					count += count_valid + sz;
					for (j = 0;j < sz;j++) {
						value = behaviorlist[t][b][dim].values[j];
						value2count[value] += 1;
					}
				}
			}
		}
	}
	if (count_pre > count) {
		return -1;
	}
	for (iter_map = value2count.begin();iter_map != value2count.end();iter_map++) {
		if (iter_map->second >= bounds[dim+1]) {
			valueset.insert(iter_map->first);
		}
	}
	sz = valueset.size();
	if (sz == 0 || sz == GetNumBehavior(tartan)) {
		tartan.dimlist[dim] = false;
		tartan.valuesetlist[dim].clear();
		return -1;
	}
	tartan.dimlist[dim] = true;
	tartan.valuesetlist[dim] = valueset;
	return 0;
}

int AddFirstSlice(struct Tartan &tartan) {
	int sz_t, t, bound;
	sz_t = tartan.tlist.size();
	t = tartan.tlist[0];
	if (t == 0) {
		return -1;
	}
	bound = (int) (tartan.blist[0].size() * _decay);
	t--;

	int b, i, j, value;
	bool isValidBehavior, isValidCell;
	set<int>::iterator iter_set;
	vector<int> bs;
	for (b = 0;b < sz_behaviors[t];b++) {
		isValidBehavior = true;
		for (i = 0;i < num_dims;i++) {
			if (tartan.dimlist[i]) {
				if (behaviorlist[t][b][i].sz == 0) {
					isValidBehavior = false;
					break;
				}
				isValidCell = false;
				for (j = 0;j < behaviorlist[t][b][i].sz;j++) {
					value = behaviorlist[t][b][i].values[j];
					iter_set = tartan.valuesetlist[i].find(value);
					if (iter_set != tartan.valuesetlist[i].end()) {
						isValidCell = true;
						break;
					}
				}
				if (!isValidCell) {
					isValidBehavior = false;
					break;
				}
			}
		}
		if (isValidBehavior) {
			bs.push_back(b);
		}
	}
	if (bs.size() < bounds[0]) {
		return -1;
	}
	tartan.tlist.insert(tartan.tlist.begin(), t);
	tartan.blist.insert(tartan.blist.begin(), bs);
	return 0;
}

int AddLastSlice(struct Tartan &tartan) {
	int sz_t, t, bound;
	sz_t = tartan.tlist.size();
	t = tartan.tlist[sz_t-1];
	if (t == num_slices-1) {
		return -1;
	}
	bound = (int) (tartan.blist[sz_t-1].size() * _decay);
	t++;

	int b, i, j, value;
	bool isValidBehavior, isValidCell;
	set<int>::iterator iter_set;
	vector<int> bs;
	for (b = 0;b < sz_behaviors[t];b++) {
		isValidBehavior = true;
		for (i = 0;i < num_dims;i++) {
			if (tartan.dimlist[i]) {
				if (behaviorlist[t][b][i].sz == 0) {
					isValidBehavior = false;
					break;
				}
				isValidCell = false;
				for (j = 0;j < behaviorlist[t][b][i].sz;j++) {
					value = behaviorlist[t][b][i].values[j];
					iter_set = tartan.valuesetlist[i].find(value);
					if (iter_set != tartan.valuesetlist[i].end()) {
						isValidCell = true;
						break;
					}
				}
				if (!isValidCell) {
					isValidBehavior = false;
					break;
				}
			}
		}
		if (isValidBehavior) {
			bs.push_back(b);
		}
	}
	if (bs.size() < bounds[0]) {
		return -1;
	}
	tartan.tlist.push_back(t);
	tartan.blist.push_back(bs);
	return 0;
}

int RemoveFirstSlice(struct Tartan &tartan) {
	int sz_t = tartan.tlist.size();
	if (sz_t <= 1 || tartan.blist[0].size() >= bounds[0]) {
		return -1;
	}
	tartan.tlist.erase(tartan.tlist.begin());
	tartan.blist.erase(tartan.blist.begin());
	return 0;
}

int RemoveLastSlice(struct Tartan &tartan) {
	int sz_t = tartan.tlist.size();
	if (sz_t <= 1 || tartan.blist[sz_t-1].size() >= bounds[0]) {
		return -1;
	}
	tartan.tlist.pop_back();
	tartan.blist.pop_back();
	return 0;
}

void *CatchTartanThread(void *id) {
	unsigned long long thread_id = (long long)id;
	unsigned long long randseed = (long long)id+1;
	int seed_id, it, i;
	int num_behavior_prev, num_behavior, count_prev;
	clock_t start, finish;

	for (seed_id = 0;seed_id < num_seeds;seed_id++) {
		cout << "Thread " << thread_id << " Seed " << seed_id << endl;

		struct Tartan tartan;
		GenerateSeed(tartan, randseed);
//		PrintTartan(tartan, "iter 0 SEED");

		num_behavior_prev = GetNumBehavior(tartan);
		count_prev = 0;
		start = clock();
		for (it = 0;it < max_iter;it++) {
			// 0. update a dimension i
			for (i = 0;i < num_dims;i++) {
				UpdateBehavior(tartan);
				UpdateDimension(tartan, i);
//				PrintTartan(tartan, "iter "+itos(it+1)+" UBD("+itos(i)+")");
			}
			// 1. check the previous slice (add?)
			AddFirstSlice(tartan);
//			PrintTartan(tartan, "iter "+itos(it+1)+" AFS");
			// 2. check the next slice (add?)
			AddLastSlice(tartan);
//			PrintTartan(tartan, "iter "+itos(it+1)+" ALS");
			// 3. check the first slice (remove?)
			RemoveFirstSlice(tartan);
//			PrintTartan(tartan, "iter "+itos(it+1)+" RFS");
			// 4. check the last slice (remove?)
			RemoveLastSlice(tartan);
//			PrintTartan(tartan, "iter "+itos(it+1)+" RLS");
			num_behavior = GetNumBehavior(tartan);
			if (num_behavior == num_behavior_prev) {
				count_prev++;
			} else {
				num_behavior_prev = num_behavior;
				count_prev = 0;
			}
			if (count_prev == 3) break;
		}
		finish = clock();
		double timecost = (double)(finish - start) / CLOCKS_PER_SEC;

		OutputTartan(tartan, timecost, it, (int) thread_id, seed_id);
	}

	pthread_exit(NULL);
}

int main(int argc, char **argv) {
	int i;
	clock_t start, finish;

	// t timestr, 0 username, 1 location
	// 2 phrases, 3 hashtags, 4 urls
	// 5 rtnames, 6 atnames
	stringstream ss; ss << argv[1];
	workdir = ss.str();
	ss.str(""); ss << argv[2];
	datafile = ss.str();
	ss.str(""); ss << argv[3];
	valuefile = ss.str();
	ss.str(""); ss << argv[4];
	outfile = ss.str();

	num_threads = stoi(argv[5]);
	num_seeds = stoi(argv[6]);
	sz_seed = stoi(argv[7]);
	max_iter = 100;

	start = clock();
	LoadData();
	finish = clock();
	double timecost = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "Timecost " << "LoadData: " << timecost << endl;

	pthread_t *pt = (pthread_t *)malloc(num_threads * sizeof(pthread_t));		
	for (i = 0;i < num_threads;i++) pthread_create(&pt[i], NULL, CatchTartanThread, (void *) i);
	for (i = 0;i < num_threads;i++) pthread_join(pt[i], NULL);

	return 0;
}


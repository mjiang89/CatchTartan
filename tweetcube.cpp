/*
TweetCube: Automatic Stream Cube Summarization for Multidimensional Social Media Mining

Contact Author: Meng Jiang, UIUC, mjiang89@illinois.edu, mjiang89@gmail.com
*/

// Format of the data file:
//
// The data file contains serveral lines, each line represents a cell and its weight in the cube.
// More specifically, each line has the following format "<d1> ... <dn> <w>", meaning a cell at <di> value on the i-th dimension with weight as <w>.
// <d1> ... <dn> and <w> are seperated by ' ' or '\t' (blank or tab)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <gsl/gsl_rng.h>

#include <map>
#include <set>
#include <vector>

#define MAX_STRING 256

typedef double real;
typedef long long bits;

struct Shape {
	int *dimsize;
	int num_entries;
	int sum_weights;
};

struct Block {
	std::set<int> *valueset;
	std::set<int> *entryset;
	struct Shape *block_shape;
};

char work_dir[MAX_STRING], input_file[MAX_STRING];
int dim = 3, metric = 2, depth = 5, seedsize = 100, num_seeds = 5, num_threads = 2, max_iter = 50;
struct Shape *data_shape;
int *entries;
std::set<int> *value2entry;
std::vector< std::pair<real, struct Block *> > *best_blocks;
bits dl_data, v_data;
FILE *fresult;

// Fast generate a random integer in [0, table_size)
int Rand(unsigned long long &randseed, int table_size) {
	randseed = randseed * 25214903917 + 11;
	return (randseed >> 16) % table_size;
}

// Log*(x) = Log(x) + Log(Log(x)) + Log(Log(Log(x))) ...
int LogStar(int x) {
	if (x <= 1) return 0;
	if (x == 2) return 1;
	int i = 0;
	while (x) {
		i++;
		x = (x >> 1);
	}
	return i + LogStar(i);
}

// Description length of [ns_1]x...x[ns_d] data with [m] 1s
void LenData() {
	dl_data = 0;
	v_data = 1;
	bits sequence;
	int i;
	real rho, _rho, h;
	std::map<int, int>::iterator it_map;
	dl_data += LogStar(data_shape->sum_weights); // for number of ones
	for (i = 0; i < dim; i++) {
		dl_data += LogStar(data_shape->dimsize[i]); // for scales
		v_data *= data_shape->dimsize[i];
	}
	sequence = v_data + data_shape->sum_weights;
	rho = 1.0 * data_shape->sum_weights / sequence;
	_rho = 1.0 - rho;
	h = 0;
	if (rho > 0) h -= rho * log2(rho);
	if (_rho > 0) h -= _rho * log2(_rho);
	dl_data += sequence * h; // for sequence
}

// Description length of a [ns_block_1]x...x[ns_block_d] block with [m_block] 1s and the rest part
// in [ns_data_1]x...x[ns_data_d] data with [m_data] 1s
bits LenBlock(struct Shape *block_shape) {
	bits ret = 0, v = 1, v_rest, sequence, sequence_rest;
	int i, sum_weights_rest;
	real rho, _rho, h;
	sum_weights_rest = data_shape->sum_weights - block_shape->sum_weights;
	// sequence block
	ret += LogStar(block_shape->sum_weights); // for number of ones
	for (i = 0; i < dim; i++) {
		ret += LogStar(block_shape->dimsize[i]); // for scales
		v *= block_shape->dimsize[i];
	}
	v_rest = v_data - v;
	sequence = v + block_shape->sum_weights;
	rho = 1.0 * block_shape->sum_weights / sequence;
	_rho = 1.0 - rho;
	h = 0;
	if (rho > 0) h -= rho * log2(rho);
	if (_rho > 0) h -= _rho * log2(_rho);
	ret += sequence * h; // for sequence
	// sequence rest
	ret += LogStar(sum_weights_rest); // for number of ones
	for (i = 0; i < dim; i++) {
		ret += LogStar(data_shape->dimsize[i]); // for scales
	}
	sequence_rest = v_rest + sum_weights_rest;
	rho = 1.0 * sum_weights_rest / sequence_rest;
	_rho = 1.0 - rho;
	h = 0;
	if (rho > 0) h -= rho * log2(rho);
	if (_rho > 0) h -= _rho * log2(_rho);
	ret += sequence_rest * h; // for sequence
	return ret;
}

// Description length of the partition as a [ns_block_1]x...x[ns_block_d] block
// in [ns_data_1]x...x[ns_data_d] data
bits LenPartition(struct Shape *block_shape) {
	bits ret = 0, sequence;
	int i;
	real ratio, ratio_rest, h;
	for (i = 0; i < dim; i++) {
		ret += LogStar(data_shape->dimsize[i]); // for data scales
		// sequence
		ret += LogStar(block_shape->dimsize[i]); // for number of 1s
		sequence = data_shape->dimsize[i];
		ratio = 1.0 * block_shape->dimsize[i] / sequence;
		ratio_rest = 1.0 - ratio;
		h = 0;
		if (ratio > 0) h -= ratio * log2(ratio);
		if (ratio_rest > 0) h -= ratio_rest * log2(ratio_rest);
		ret += sequence * h; // for sequence
	}
	return ret;
}

// Suspiciousness with ER-Poisson probability principle [ICDM'15]
real Suspiciousness(struct Shape *block_shape) {
	int i;
	real temp1 = 1, temp2 = 0, p;
	for (i = 0; i < dim; i++) {
		p = 1.0 * block_shape->dimsize[i] / data_shape->dimsize[i];
		temp1 *= p;
		temp2 += log(p);
	}
	temp1 *= data_shape->sum_weights;
	temp2 *= block_shape->sum_weights;
	return log(block_shape->sum_weights * (log(1.0 * block_shape->sum_weights / data_shape->sum_weights) - 1.0) + temp1 - temp2);
}

// Read data from file with [postfix]
int ReadData(char *postfix) {
	int i, j, k, value, weight, offset, sum_dimsize;
	char data_file[MAX_STRING];
	strcpy(data_file, work_dir);
	strcat(data_file, "/");
	strcat(data_file, input_file);
	strcat(data_file, postfix);
	printf("data_file: %s\n",data_file);

	data_shape = new struct Shape();
	data_shape->dimsize = (int *)malloc(dim * sizeof(int));
	data_shape->num_entries = 0;
	data_shape->sum_weights = 0;

	FILE *fin;
	char name[MAX_STRING], line[MAX_STRING];
	fin = fopen(data_file, "rb");
	if (fin == NULL) {
		printf("ERROR: input file not found!\n");
		return -1;
	}
	while (fgets(line, sizeof(line), fin)) data_shape->num_entries++;
	if (data_shape->num_entries == 0) {
		printf("ERROR: input file is empty!\n");
		return -1;
	}
	fclose(fin);
	// Test [num_entries]
	printf("  Number of nonzero entries: %d\n", data_shape->num_entries);

	for (i = 0; i < dim; i++) data_shape->dimsize[i] = 0;
	sum_dimsize = 0;
	fin = fopen(data_file, "rb");
	for (k = 0; k < data_shape->num_entries; k++) {
		for (i = 0; i < dim; i++) {
			fscanf(fin, "%s", name);
			value = atoi(name);
			if (value > data_shape->dimsize[i]) data_shape->dimsize[i] = value;
		}
		fscanf(fin, "%s", name);
		weight = atoi(name);
		data_shape->sum_weights += weight;
	}
	fclose(fin);
	for (i = 0; i < dim; i++) {
		data_shape->dimsize[i]++;
		sum_dimsize += data_shape->dimsize[i];
	}
	// Test [dimsize] [sum_weights]
	for (i = 0; i < dim; i++) printf("  Size of dimension [%d]: %d\n", i, data_shape->dimsize[i]);
	printf("  Sum of weights: %d\n", data_shape->sum_weights);

	entries = (int *)malloc(data_shape->num_entries * (dim + 1) * sizeof(int));
	value2entry = new std::set<int>[sum_dimsize];
	fin = fopen(data_file, "rb");
	for (k = 0; k < data_shape->num_entries; k++) {
		for (i = 0; i < dim; i++) {
			fscanf(fin, "%s", name);
			value = atoi(name);
			entries[k * (dim + 1) + i] = value;
			offset = 0;
			for (j = 0; j < i; j++) offset += data_shape->dimsize[j];
			value2entry[offset + value].insert(k);
		}
		fscanf(fin, "%s", name);
		weight = atoi(name);
		entries[k * (dim + 1) + dim] =  weight;
	}
	fclose(fin);

	best_blocks = new std::vector< std::pair<real, struct Block *> >();
	return 0;
}

// Generate a seed block with [randseed] from thread
void GenerateSeedBlock(struct Block *block, unsigned long long &randseed) {
	int t, i, j, k, entry_id, value, weight, offset;
	std::set<int> temp_set;
	std::set<int>::iterator it_set;
	for (t = 0; t < seedsize; t++) {
		entry_id = Rand(randseed, data_shape->num_entries);
		for (i = 0; i < dim; i++) {
			value = entries[entry_id * (dim + 1) + i];
			block->valueset[i].insert(value);
			offset = 0;
			for (j = 0; j < i; j++) offset += data_shape->dimsize[j];
			for (it_set = value2entry[offset + value].begin(); it_set != value2entry[offset + value].end(); it_set++) {
				k = *it_set;
				block->entryset[i].insert(k);
			}

		}
	}
	std::set<int> block_entryset = block->entryset[0];
	for (i = 1; i < dim; i++) {
		temp_set.clear();
		std::set_intersection(block_entryset.begin(), block_entryset.end(),
			block->entryset[i].begin(), block->entryset[i].end(),
			std::inserter(temp_set, temp_set.begin()));
		block_entryset = temp_set;
	}
	block->block_shape->num_entries = block_entryset.size();
	for (i = 0; i < dim; i++) block->block_shape->dimsize[i] = block->valueset[i].size();
	for (it_set = block_entryset.begin(); it_set != block_entryset.end(); it_set++) {
		k = *it_set;
		weight = entries[k * (dim + 1) + dim];
		block->block_shape->sum_weights += weight;
	}
}

// Compute saving ratio of a given [block] in the data: (dl_data - dl_block - dl_partition) / dl_data
real ComputeSaving(struct Block *block) {
	real ret = 0;
	bits dl_block, dl_partition;
	if (metric == 1) {
		ret = Suspiciousness(block->block_shape);
	}
	if (metric == 2) {
		dl_block = LenBlock(block->block_shape);
		dl_partition = LenPartition(block->block_shape);
		ret = 1.0 * (dl_data - dl_block - dl_partition) / dl_data;
	}
	return ret;
}

// Comparer function for sorting "map" (vector<pair>) by values 
bool value_comparer(const std::pair<int, int> &a, const std::pair<int, int> &b) {
	return a.second > b.second;
}

// Update a [block] at dimension [dim_curr]
real UpdateBlockDimension(struct Block *block, int dim_curr) {
	int i, j, k, value, weight, offset;
	real saving_curr, saving_max;
	std::set<int> block_entryset;
	std::set<int> block_entryset_other;
	std::set<int> temp_set;
	std::set<int> entryset_prev, block_entryset_prev;
	std::map<int, int> value2sumweight;
	std::vector< std::pair<int, int> > value_sumweight;
	std::set<int>::iterator it_set, _it_set;
	std::map<int, int>::iterator it_map;
	std::vector< std::pair<int, int> >::iterator it_vector;

	j = 0;
	for (i = 0; i < dim; i++) {
		if (i != dim_curr) {
			if (j == 0) {
				block_entryset_other = block->entryset[i];
			} else {
				temp_set.clear();

				std::set_intersection(block_entryset_other.begin(), block_entryset_other.end(),
					block->entryset[i].begin(), block->entryset[i].end(),
				std::inserter(temp_set, temp_set.begin()));
				/* 2nd option
				if (block_entryset_other.size() < block->entryset[i].size()) {
					for (it_set = block_entryset_other.begin(); it_set != block_entryset_other.end(); it_set++) {
						k = *it_set;
						_it_set = block->entryset[i].find(k);
						if (_it_set != block->entryset[i].end()) {
							temp_set.insert(k);
						}
					}
				} else {
					for (it_set = block->entryset[i].begin(); it_set != block->entryset[i].end(); it_set++) {
						k = *it_set;
						_it_set = block_entryset_other.find(k);
						if (_it_set != block_entryset_other.end()) {
							temp_set.insert(k);
						}
					}
				}
				*/

				block_entryset_other = temp_set;
			}
			j = 1;
		}
	}
	for (it_set = block_entryset_other.begin(); it_set != block_entryset_other.end(); it_set++) {
		k = *it_set;
		value = entries[k * (dim + 1) + dim_curr];
		weight = entries[k * (dim + 1) + dim];
		value2sumweight[value] += weight;
	}
	for (it_map = value2sumweight.begin(); it_map != value2sumweight.end(); it_map++) {
		value_sumweight.push_back(*it_map);
	}
	std::sort(value_sumweight.begin(), value_sumweight.end(), value_comparer);

	// Initial by the larget-weighted value on the dimension
	it_vector = value_sumweight.begin();
	value = it_vector->first;
	block->valueset[dim_curr].clear();
	block->valueset[dim_curr].insert(value);
	offset = 0;
	for (j = 0; j < dim_curr; j++) offset += data_shape->dimsize[j];
	block->entryset[dim_curr] = value2entry[offset + value];
	block->block_shape->dimsize[dim_curr] = 1;
	block->block_shape->num_entries = 0;
	block->block_shape->sum_weights = 0;
	for (it_set = block->entryset[dim_curr].begin(); it_set != block->entryset[dim_curr].end(); it_set++) {
		k = *it_set;
		_it_set = block_entryset_other.find(k);
		if (_it_set != block_entryset_other.end()) {
			block_entryset.insert(k);
			weight = entries[k * (dim + 1) + dim];
			block->block_shape->num_entries += 1;
			block->block_shape->sum_weights += weight;
		}
	}
	saving_max = ComputeSaving(block);

	// Find the best set of values on the dimension
	it_vector++;
	while (it_vector != value_sumweight.end()) {
		value = it_vector->first;
		block->valueset[dim_curr].insert(value);
		block->block_shape->dimsize[dim_curr]++;
		offset = 0;
		for (j = 0; j < dim_curr; j++) offset += data_shape->dimsize[j];
		for (it_set = value2entry[offset + value].begin(); it_set != value2entry[offset + value].end(); it_set++) {
			k = *it_set;
			block->entryset[dim_curr].insert(k);
			entryset_prev.insert(k);
			_it_set = block_entryset_other.find(k);
			if (_it_set != block_entryset_other.end()) {
				block_entryset.insert(k);
				block_entryset_prev.insert(k);
				weight = entries[k * (dim + 1) + dim];
				block->block_shape->num_entries += 1;
				block->block_shape->sum_weights += weight;
			}
		}
		saving_curr = ComputeSaving(block);

		if (saving_curr > saving_max) {
			saving_max = saving_curr;
		} else {
			it_set = block->valueset[dim_curr].find(value);
			block->valueset[dim_curr].erase(it_set);
			for (it_set = entryset_prev.begin(); it_set != entryset_prev.end(); it_set++) {
				k = *it_set;
				_it_set = block->entryset[dim_curr].find(k);
				block->entryset[dim_curr].erase(_it_set);
			}
			block->block_shape->dimsize[dim_curr]--;
			for (it_set = block_entryset_prev.begin(); it_set != block_entryset_prev.end(); it_set++) {
				k = *it_set;
				_it_set = block_entryset.find(k);
				block_entryset.erase(_it_set);
				weight = entries[k * (dim + 1) + dim];
				block->block_shape->num_entries -= 1;
				block->block_shape->sum_weights -= weight;
			}
			break;
		}
		entryset_prev.clear();
		block_entryset_prev.clear();
		it_vector++;
	}
	// Test "after adjust one dimension"
	printf("\t  Operate dimension [%d]: selecting from %lu values\n", dim_curr, value2sumweight.size());
	printf("\t   Number of nonzero entries: %d\n", block->block_shape->num_entries);
	for (i = 0; i < dim; i++) printf("\t    Size of dimension [%d]: %d / %d\n", i, block->block_shape->dimsize[i], data_shape->dimsize[i]);
	printf("\t    Sum of weights: %d\n", block->block_shape->sum_weights);

	return saving_max;
}

// Thread that finds dense blocks with [num_seeds] seeds at the current depth
void *FindBlockThread(void *id) {
	int t, i, d, d_prev, dim_prev, dim_curr, best_t = -1, best_iter = -1;
	unsigned long long thread_id = (long long)id;
	unsigned long long randseed = (long long)id;
	real saving, saving_prev = -1.0, best_saving = -1.0;
	struct Block *best_block;
	clock_t start, finish;
	clock_t start_seed, finish_seed;

	start = clock();
	for (t = 0; t < num_seeds; t++) {
		struct Block *block = new Block();
		block->valueset = new std::set<int>[dim];
		block->entryset = new std::set<int>[dim];
		block->block_shape = new struct Shape();
		block->block_shape->dimsize = (int *)malloc(dim * sizeof(int));
		block->block_shape->num_entries = 0;
		block->block_shape->sum_weights = 0;

		GenerateSeedBlock(block, randseed);
		// Test "GenerateSeedBlock"
		printf("\tThread %llu: seed block #%d\n", thread_id, t);
		printf("\t  Number of nonzero entries: %d\n", block->block_shape->num_entries);
		for (i = 0; i < dim; i++) printf("\t    Size of dimension [%d]: %d / %d\n", i, block->block_shape->dimsize[i], data_shape->dimsize[i]);
		printf("\t  Sum of weights: %d\n", block->block_shape->sum_weights);

		dim_prev = -1;
		d_prev = -1;
		for (d = 0; d < max_iter; d++) {
			if (d - d_prev > dim) break;

			start_seed = clock();
			while (1) {
				dim_curr = Rand(randseed, dim);
				if (dim_curr != dim_prev) break;
			}

			// Test "UpdateBlockDimension" (begin)
			printf("\t-- Thread %llu, Seed #%d, Iter #%d begins: dimension [%d] --\n", thread_id, t, d, dim_curr);

			dim_prev = dim_curr;
			saving = UpdateBlockDimension(block, dim_curr);
			if (saving != saving_prev) {
				saving_prev = saving;
				d_prev = d;
			}
			finish_seed = clock();
			
			// Test "UpdateBlockDimension" (end)
			printf("\t-- Thread %llu, Seed #%d, Iter #%d ends: dimension [%d] => saving %.5f, time cost %lf --\n",
				thread_id, t, d, dim_curr, saving, (double)(finish_seed - start_seed) / CLOCKS_PER_SEC);
			
		}
		printf("\t-- Thread %llu, Seed #%d ends at Iter #%d --\n", thread_id, t, d);
		if (saving > best_saving) {
			best_saving = saving;
			best_block = block;
			best_t = t;
			best_iter = d;
		}
	}
	best_blocks->push_back(std::make_pair(best_saving, best_block));
	finish = clock();
	// Test "The Best Block of the thread"
	printf("\tThread %llu Best: Seed %d at Iter #%d, saving %.5f, time cost %lf\n", thread_id,
		best_t, best_iter, best_saving,
		(double)(finish - start) / CLOCKS_PER_SEC);

	pthread_exit(NULL);
}

// Output the block (1) and the rest of data (0) into files
int OutputData(char *postfix, struct Block *block) {
	int i, j, k, value;
	std::set<int>::iterator it_set;

	char data_file[MAX_STRING];
	strcpy(data_file, work_dir);
	strcat(data_file, "/");
	strcat(data_file, input_file);
	strcat(data_file, postfix);
	char data_file0[MAX_STRING];
	strcpy(data_file0, data_file);
	strcat(data_file0, "0");
	char data_file1[MAX_STRING];
	strcpy(data_file1, data_file);
	strcat(data_file1, "1");

	FILE *fin, *fout0, *fout1;
	char name[MAX_STRING], line[MAX_STRING];
	fin = fopen(data_file, "rb");
	if (fin == NULL) {
		printf("ERROR: input file not found!\n");
		return -1;
	}
	fout0 = fopen(data_file0, "wb");
	if (fout0 == NULL) {
		printf("ERROR: output file 0 cannot write!\n");
		return -1;
	}
	fout1 = fopen(data_file1, "wb");
	if (fout1 == NULL) {
		printf("ERROR: output file 1 cannot write!\n");
		return -1;
	}
	for (k = 0; k < data_shape->num_entries; k++) {
		j = 1;
		for (i = 0; i < dim; i++) {
			fscanf(fin, "%s", name);
			if (i == 0) {
				strcpy(line, name);
			} else {
				strcat(line, " ");
				strcat(line, name);
			}
			value = atoi(name);
			it_set = block->valueset[i].find(value);
			if (it_set == block->valueset[i].end()) j = 0;
		}
		fscanf(fin, "%s", name);
		strcat(line, " ");
		strcat(line, name);
		if (j == 0) {
			fprintf(fout0, "%s\n", line);
		} else {
			fprintf(fout1, "%s\n", line);
		}
	}
	fclose(fout1);
	fclose(fout0);
	fclose(fin);
	return 0;
}

// Tree process to find blocks
void FindBlocks(char *postfix) {
	int cur_depth = strlen(postfix);
	printf("-- Current depth: %d ('%s') --\n", cur_depth, postfix);
	if (cur_depth == depth) return;
	int i;
	clock_t start, finish;

	printf("-- (1) ReadData '%s' --\n", postfix);
	start = clock();
	if (ReadData(postfix) < 0) return;
	finish = clock();
	printf("-- (1) ReadData '%s': time cost %lf --\n", postfix, (double)(finish - start) / CLOCKS_PER_SEC);

	printf("-- (2) LenData '%s' --\n", postfix);
	start = clock();
	if (metric == 1) { // suspiciousness [ICDM'15] 
	}
	if (metric == 2) { // mdl
		LenData();
		printf("  Description length of data: %lld bits\n", dl_data);
	}
	finish = clock();
	printf("-- (2) LenData '%s': time cost %lf --\n", postfix, (double)(finish - start) / CLOCKS_PER_SEC);

	printf("-- (3) FindBlockThread '%s' --\n", postfix);
	start = clock();
	pthread_t *pt = (pthread_t *)malloc(num_threads * sizeof(pthread_t));		
	for (i = 0; i < num_threads; i++) pthread_create(&pt[i], NULL, FindBlockThread, (void *)i);
	for (i = 0; i < num_threads; i++) pthread_join(pt[i], NULL);
	std::vector< std::pair<real, struct Block *> >::iterator it_vector;
	real saving, best_saving = -1.0;
	struct Block *best_block;
	for (it_vector = best_blocks->begin(); it_vector != best_blocks->end(); it_vector++) {
		saving = it_vector->first;
		if (saving > best_saving) {
			best_saving = saving;
			best_block = it_vector->second;
		}
	}
	if (best_saving < 0) return;

	// Test "The Best Block"
	printf("The Best Block:\n");
	printf("  Number of nonzero entries: %d\n", best_block->block_shape->num_entries);
	for (i = 0; i < dim; i++) printf("  Size of dimension [%d]: %d\n", i, best_block->block_shape->dimsize[i]);
	printf("  Sum of weights: %d\n", best_block->block_shape->sum_weights);
	printf("  Best saving (suspiciousness) ratio: %.5f\n", best_saving);
	finish = clock();
	printf("-- (3) FindBlockThread '%s': time cost %lf --\n", postfix, (double)(finish - start) / CLOCKS_PER_SEC);

	if (data_shape->sum_weights == best_block->block_shape->sum_weights) return;

	printf("-- (4) OutputData '%s' --\n", postfix);
	start = clock();
	if (OutputData(postfix, best_block) < 0) return;
	finish = clock();
	printf("-- (4) OutputData: time cost %lf --\n", (double)(finish - start) / CLOCKS_PER_SEC);
	
	char _postfix[MAX_STRING];
	int num_entries_rest = data_shape->num_entries - best_block->block_shape->num_entries;
	int sum_weights_rest = data_shape->sum_weights - best_block->block_shape->sum_weights;
//	bool is_result = true;
	if (num_entries_rest >= seedsize && sum_weights_rest >= 1 * num_entries_rest) {
		// Find blocks in "-0": from the rest of data
		strcpy(_postfix, postfix);
		strcat(_postfix, "0");
		FindBlocks(_postfix);
	} else {
//		is_result = false;
	}
	if (best_block->block_shape->num_entries >= seedsize && 
		best_block->block_shape->sum_weights >= 1 * best_block->block_shape->num_entries) {
		printf("==> %d %d\n",best_block->block_shape->num_entries,seedsize);
		// Find blocks in "-1": from the block
		strcpy(_postfix, postfix);
		strcat(_postfix, "1");
		FindBlocks(_postfix);
	} else {
//		is_result = false;
	}
//	if (is_result) {
		// for result
		fprintf(fresult, "%s1", postfix);
		for (i = 0; i < dim; i++) fprintf(fresult, " %d", best_block->block_shape->dimsize[i]);
		fprintf(fresult, " %d", best_block->block_shape->num_entries);
		fprintf(fresult, " %d", best_block->block_shape->sum_weights);
		fprintf(fresult, " %.5f\n", best_saving);
//	}
}

// Main process
void RunTweetCube() {
	printf("-- Parameters --\n");
	printf("work_dir: %s\n", work_dir);
	printf("input_file: %s\n", input_file);
	printf("dim: %d\n", dim);
	printf("metric: %d\n", metric);
	printf("depth: %d\n", depth);
	printf("seedsize: %d\n", seedsize);
	printf("num_seeds: %d\n", num_seeds);
	printf("num_threads: %d\n", num_threads);
	printf("max_iter: %d\n", max_iter);

	char result_file[MAX_STRING];
	strcpy(result_file, "result_");
	strcat(result_file, input_file);
	strcat(result_file, ".txt");

	fresult = fopen(result_file, "wb");
	fprintf(fresult, "work_dir: %s\n", work_dir);
	fprintf(fresult, "input_file: %s\n", input_file);
	fprintf(fresult, "dim: %d\n", dim);
	fprintf(fresult, "metric: %d\n", metric);
	fprintf(fresult, "depth: %d\n", depth);
	fprintf(fresult, "seedsize: %d\n", seedsize);
	fprintf(fresult, "num_seeds: %d\n", num_seeds);
	fprintf(fresult, "num_threads: %d\n", num_threads);
	fprintf(fresult, "max_iter: %d\n", max_iter);
	FindBlocks("");
	fclose(fresult);
}

int ArgPos(char *str, int argc, char **argv) {
	int i;
	for (i = 1; i < argc; i++) {
		if (!strcmp(str, argv[i])) {
			if (i == argc - 1) {
				printf("Argument missing for %s\n", str);
				exit(1);
			}
			return i;
		}
	}
	return -1;
}

int main(int argc, char **argv) {
	int i;
	if (argc == 1) {
		printf("TweetCube: Automatic Stream Cube Summarization for Multidimensional Social Media Mining\n\n");
		printf("Options:\n");
		printf("Parameters:\n");
		printf("\t-workdir <dir>\n");
		printf("\t\tUse <dir> as workspace\n");
		printf("\t-input <file>\n");
		printf("\t\tUse datacube from <file> to find dense blocks\n");
		printf("\t-dim <int>\n");
		printf("\t\tSet the number of dimensions; default is 3\n");
		printf("\t-metric <int>\n");
		printf("\t\tSet the metric type of evaluating dense blocks; 1 for suspiciousness [ICDM'15], 2 for mdl; default is 2\n");
		printf("\t-depth <int>\n");
		printf("\t\tSet the depth of searching blocks; default is 5\n");
		printf("\t-seedsize <int>\n");
		printf("\t\tSet the number of nonzero entries in a seed; default is 100\n");
		printf("\t-seeds <int>\n");
		printf("\t\tSet the number of seeds per thread as <int>; default is 5\n");
		printf("\t-threads <int>\n");
		printf("\t\tUse <int> threads; default is 2\n");
		printf("\t-iter <int>\n");
		printf("\t\tSet the maximum number of iterations; default is 50\n");
		printf("\nExamples:\n");
		printf("./tweetcube -workdir data -input localevent -dim 3 -metric 2 -depth 5 -seedsize 100 -seeds 5 -threads 2 -iter 50\n");
		return 0;
	}
	if ((i = ArgPos((char *)"-workdir", argc, argv)) > 0) strcpy(work_dir, argv[i + 1]);
	if ((i = ArgPos((char *)"-input", argc, argv)) > 0) strcpy(input_file, argv[i + 1]);
	if ((i = ArgPos((char *)"-dim", argc, argv)) > 0) dim = atoi(argv[i + 1]);
	if ((i = ArgPos((char *)"-metric", argc, argv)) > 0) metric = atoi(argv[i + 1]);
	if ((i = ArgPos((char *)"-depth", argc, argv)) > 0) depth = atoi(argv[i + 1]);
	if ((i = ArgPos((char *)"-seedsize", argc, argv)) > 0) seedsize = atoi(argv[i + 1]);
	if ((i = ArgPos((char *)"-seeds", argc, argv)) > 0) num_seeds = atoi(argv[i + 1]);
	if ((i = ArgPos((char *)"-threads", argc, argv)) > 0) num_threads = atoi(argv[i + 1]);
	if ((i = ArgPos((char *)"-iter", argc, argv)) > 0) max_iter = atoi(argv[i + 1]);

	RunTweetCube();
	return 0;
}



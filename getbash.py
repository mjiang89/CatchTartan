import itertools

fw = open('compile.sh','w')
fw.write('rm -f catchtartan\n')
fw.write('g++ -std=c++0x -lm -pthread -O3 -march=native -Wall -funroll-loops -ffast-math -Wno-unused-result catchtartan.cpp -o catchtartan -lgsl -lm -lgslcblas\n')
fw.write('g++ -std=c++0x -lm -pthread -O3 -march=native -Wall -funroll-loops -ffast-math -Wno-unused-result catchtartan-crj.cpp -o catchtartan-crj -lgsl -lm -lgslcblas\n')
fw.close()

workdir = '/usr0/home/mjiang89/'
num_threads = 10
num_seeds = 30
sz_seed = 300
for dataset in ['LA14','NYC14','SPB13','GRM13']:
	fw = open('run_'+dataset+'.sh','w')
	fw.write('mkdir output_'+dataset+'\n')
	s = './catchtartan '+workdir+' data/data_'+dataset+'.txt value/value_'+dataset+'_ output_'+dataset+'/output_'+dataset+'_ '+str(num_threads)+' '+str(num_seeds)+' '+str(sz_seed)
	fw.write(s+'\n')
	fw.close()

workdir = '/Users/meng/Desktop/tartan/'
num_threads = 2
num_seeds = 5
sz_seed = 300

for dataset in ['CRJ']:
	fw = open('run_'+dataset+'.sh','w')
	fw.write('mkdir output_'+dataset+'\n')
	s = './catchtartan-crj '+workdir+' data/data_'+dataset+'.txt value/value_'+dataset+'_ output_'+dataset+'/output_'+dataset+'_ '+str(num_threads)+' '+str(num_seeds)+' '+str(sz_seed)
	fw.write(s+'\n')
	fw.close()



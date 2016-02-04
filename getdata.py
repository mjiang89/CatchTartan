import sys

# 0: timestr
# 1: username
# 2: location
# 3: list of phrases
# 4: list of hashtags
# 5: list of urls
# 6: list of rtnames
# 7: list of atnames

# 8: timestamp
# 9: tweet

def GetTweet(dataset):
	stopwordset = set()
	fr = open('tweet/stopword.txt','rb')
	for line in fr:
		word = line.strip('\r\n')
		stopwordset.add(word)
	fr.close()
	fw = open('tweet/tweet_'+dataset+'.txt','w')
	fr = open('tweet/tweet0_'+dataset+'.txt','rb')
	for line in fr:
		arr = line.strip('\r\n').split('|')
		s = ''
		for i in range(0,len(arr)):
			if i == 3 and len(arr[i]) > 0:
				item = ''
				arr0 = arr[i].split(',')
				for x in arr0:
					phrase = x[1:len(x)-1]
					isnotstop = False
					arr1 = phrase.split(' ')
					for word in arr1:
						if len(word) > 0 and not word in stopwordset:
							isnotstop = True
							break
					if isnotstop:
						item += ','+x
				if len(item) > 0: item = item[1:]
				s += '|'+item
			else:
				s += '|'+arr[i]
		s = s[1:]
		fw.write(s+'\n')
	fr.close()
	fw.close()

def GetValue(dataset):
	entry2counts = [{} for i in range(0,8)]
	fr = open('tweet/tweet_'+dataset+'.txt','rb')
	for line in fr:
		arr = line.strip('\r\n').split('|')
		for i in range(0,3):
			entry = arr[i]
			if len(entry) == 0: continue
			if not entry in entry2counts[i]:
				entry2counts[i][entry] = 0
			entry2counts[i][entry] += 1
		for i in range(3,8):
			for entry in arr[i].split(','):
				if len(entry) == 0: continue
				if not entry in entry2counts[i]:
					entry2counts[i][entry] = 0
				entry2counts[i][entry] += 1
	fr.close()
	for i in range(0,8):
		n = len(entry2counts[i])
		fw = open('value/value_'+dataset+'_'+str(i)+'.txt','w')
		fw.write('#'+str(n)+'\n')
		if i == 0:
			j = -1
			for e in sorted(entry2counts[i].items(),key=lambda x:x[0]):
				j += 1
				fw.write(str(j)+'|'+e[0]+'|'+str(e[1])+'\n')
		else:
			j = -1
			for e in sorted(entry2counts[i].items(),key=lambda x:-x[1]):
				j += 1
				fw.write(str(j)+'|'+e[0]+'|'+str(e[1])+'\n')
		fw.close()

def GetData(dataset):
	fw = open('data/data_'+dataset+'.txt','w')
	s = ''
	entry2values = [{} for i in range(0,8)]
	for i in range(0,8):
		fr = open('value/value_'+dataset+'_'+str(i)+'.txt','rb')
		line = fr.readline()
		s += '|'+line.strip('\r\n')
		for line in fr:
			arr = line.strip('\r\n').split('|')
			entry,value = arr[1],int(arr[0])
			entry2values[i][entry] = value
		fr.close()
	fw.write(s[1:]+'\n')
	fr = open('tweet/tweet_'+dataset+'.txt','rb')
	for line in fr:
		arr = line.strip('\r\n').split('|')
		s = ''
		for i in range(0,3):
			value = ''
			entry = arr[i]
			if not len(entry) == 0 and entry in entry2values[i]:
				value = str(entry2values[i][entry])
			s += '|'+value
		for i in range(3,8):
			values = ''
			for entry in arr[i].split(','):
				value = ''
				if not len(entry) == 0 and entry in entry2values[i]:
					value = str(entry2values[i][entry])
				values += ','+value
			if len(values) > 0: values = values[1:]
			s += '|'+values
		if len(s) > 0: s = s[1:]
		fw.write(s+'\n')
		fw.write(arr[8]+'|'+arr[1]+'|'+arr[2]+'|'+arr[9]+'\n')
	fr.close()
	fw.close()

def Check():
	fr = open('data/data_NYC14.txt','rb')
	fr.readline()
	line = fr.readline()
	while len(line) > 0:
		arr = line.strip('\r\n').split('|')
		arr0 = arr[3].split(',')
		valueset = set()
		isvalid = False
		for x in arr0:
			if len(x) > 0:
				value = int(x)
				if value in valueset:
					print line.strip('\r\n')
					line = fr.readline()
					print line.strip('\r\n')
					isvalid = True
					break
				else:
					valueset.add(value)
		if isvalid: break
		fr.readline()
		line = fr.readline()
	fr.close()

def GetNumDim2Count():
	fw = open('numdim2count.csv','w')
	for dataset in ['NYC14','LA14','SPB13','GRM13']:
		numdim2count = {}
		numvalue2count = {}
		totalcount,totalcount_dim,totalcount_value = 0,0,0
		fr = open('data/data_'+dataset+'.txt','rb')
		fr.readline()
		line = fr.readline()
		while len(line) > 0:
			arr = line.strip('\r\n').split('|')
			numdim,numvalue = 0,0
			for i in range(1,len(arr)):
				if len(arr[i]) > 0:
					numdim += 1
					arr0 = arr[i].split(',')
					numvalue += len(arr0)
			if not numdim in numdim2count:
				numdim2count[numdim] = 0
			numdim2count[numdim] += 1
			if not numvalue in numvalue2count:
				numvalue2count[numvalue] = 0
			numvalue2count[numvalue] += 1
			fr.readline()
			line = fr.readline()
		fr.close()
		for numdim in numdim2count:
			totalcount += numdim2count[numdim]
			totalcount_dim += numdim*numdim2count[numdim]
		for numvalue in numvalue2count:
			totalcount_value += numvalue*numvalue2count[numvalue]
		fw.write('#'+dataset+','+str(totalcount)+','+str(totalcount_dim)+','+str(totalcount_value)+','+str(1.0*totalcount_dim/totalcount)+','+str(1.0*totalcount_value/totalcount)+'\n')
		for e in sorted(numdim2count.items(),key=lambda x:x[0]):
			fw.write(dataset+','+str(e[0])+','+str(e[1])+','+str(1.0*e[1]/totalcount)+'\n')
		for e in sorted(numvalue2count.items(),key=lambda x:x[0]):
			fw.write(dataset+','+str(e[0])+','+str(e[1])+','+str(1.0*e[1]/totalcount)+'\n')
	fw.close()


if __name__ == '__main__':
	'''
	for dataset in ['NYC14','LA14','SPB13','GRM13']:
		GetTweet(dataset)
		GetValue(dataset)
		GetData(dataset)
	Check()
	GetNumDim2Count()
	'''

	wordset = set(['predict','will win','will lose', \
		'7-3','3-7','21-6','6-21', \
		'half','destiny','beyonce', \
		'28-23','23-28','34-31','31-34', \
		'champion','congratulation','game over'])
	fw = open('tweet/tweet_CRJ.txt','w')
	fr = open('tweet/tweet_SPB13.txt','rb')
	for line in fr:
		text = line.strip('\r\n')
		isValid = False
		for word in wordset:
			if word in text:
				isValid = True
				break
		if isValid:
			fw.write(text+'\n')
	fr.close()
	fw.close()
	GetValue('CRJ')
	GetData('CRJ')


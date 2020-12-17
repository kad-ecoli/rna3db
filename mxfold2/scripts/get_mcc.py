#! /usr/bin/python3
import sys
import pandas as pd

def main():
	df = pd.read_csv(sys.argv[1])
	desc = df.describe()
	mcc = desc.loc['mean']['mcc']
	print(mcc)

if __name__ == '__main__':
	main()

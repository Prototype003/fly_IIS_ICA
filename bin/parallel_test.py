import multiprocessing

def parallel_function(fly):
	for i in range(0, 10):
		print(fly)

def main():
	flies = [1, 2, 3, 4]
	for fly in flies:
		process = multiprocessing.Process(target=parallel_function, args=(fly,))
		process.start()
		process.join()
	print('done')

if __name__ == "__main__":
	main()

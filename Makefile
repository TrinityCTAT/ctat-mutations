all:
	./download_cromwell.sh
	cd plugins/bamsifter/ && make

TARGETS = libs


# Make and install all the libraries:
libs:	
	$(MAKE) -C alps_src
	$(MAKE) -C alps_src instlocal
	


# Cleans all the library source directors and empties ./include and ./lib:
clean:
	$(MAKE) -C alps_src clean
	rm -rf include/*
	rm -f lib/*



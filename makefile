CC = g++
CFLAGS = -Wall
SRCS = main.cpp readfile.cpp interpolate.cpp transform.cpp
OBJS = $(SRCS:.cpp=.o)
MAIN = main

default: $(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(MAIN)

$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) -c $(SRCS)

clean:
	$(RM) *.o $(MAIN)

install:
	@echo "Checking for essential directories"
	@test -d "data_raw" || (mkdir "data_raw"; echo "  Creating data_raw/")
	@test -d "gridlists" || (mkdir "gridlists"; echo "  Creating gridlists/")
	@test -d "output" || (mkdir "output"; echo "  Creating output/")
	@test -d "data_trimmed" || (mkdir "data_trimmed"; echo "  Creating data_trimmed/")

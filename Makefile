CXX      = c++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

test: tests/test_poisson tests/test_joukowsky
	./tests/test_poisson
	./tests/test_joukowsky

tests/test_poisson: tests/test_poisson.cpp poisson.cpp poisson.h
	$(CXX) $(CXXFLAGS) tests/test_poisson.cpp poisson.cpp -o tests/test_poisson

tests/test_joukowsky: tests/test_joukowsky.cpp joukowsky.cpp joukowsky.h
	$(CXX) $(CXXFLAGS) tests/test_joukowsky.cpp joukowsky.cpp -o tests/test_joukowsky

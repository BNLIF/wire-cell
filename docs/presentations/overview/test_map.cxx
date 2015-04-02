#include <map>
#include <iostream>
int main()
{
    typedef std::map<int,float> SparseHist; // save typing
    SparseHist sh;
    sh[2] = 20;
    sh[4] += 2; // default value springs into life
    SparseHist::iterator it, done = sh.end();
    for (it = sh.begin(); it != done; ++it) {
	std::cout << "bin #" << it->first 
		  << " content:" << it->second << std::endl;
    }
}

#include <set>
#include <iostream>
int main()
{
    typedef std::set<int> IntBag;
    IntBag ib;
    ib.insert(5);
    ib.insert(3);
    ib.insert(7);
    IntBag::iterator it, done = ib.end();
    for (it = ib.begin(); it != done; ++it) {
	std::cout << *it << std::endl;
    }
}

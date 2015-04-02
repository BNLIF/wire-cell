#include <list>
#include <iostream>
int main() {
    typedef std::list<int> MyList;
    int dat[] = { 1,2,3 };
    MyList lst(dat,dat+3);
    lst.push_front(24);  // prepend, also pop_front()
    lst.push_back(42);  // append, also pop_back()
    MyList::iterator it, done = lst.end();
    for (it = lst.begin(); it != done; ++it) {
	*it *= 100;
	std::cout << *it << std::endl;
    }
    return 0;
}

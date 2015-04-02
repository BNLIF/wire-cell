#include <vector>
int main() {
    int dat[] = {5, 10, 15};
    std::vector<int> vec(dat, dat+3);
    vec.push_back(42);  // append
    for (std::size_t ind=0; ind < vec.size(); ++ind) {
	vec[ind] *= ind+1;
	// or
	vec.at(ind) *= ind+1;
    }
    return 0;
}


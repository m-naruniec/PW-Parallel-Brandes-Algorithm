#include <fstream>
#include <map>
#include <vector>
#include <queue>
#include <stack>
#include <cassert>
#include <utility>
#include <thread>
#include <mutex>
//#include <iostream>


int THREAD_AMOUNT;
char* INPUT_NAME;
char* OUTPUT_NAME;

std::fstream file_stream;

std::vector<std::thread> thread_vector;

std::map<int, unsigned long> renaming;
std::vector<std::vector<unsigned long> > neighbours;

std::vector<double> global_BC;

unsigned long n;

unsigned long next_unprocessed = 0;

std::mutex next_mut;
std::mutex BC_mut;

bool has_neighbours(unsigned long x) {
    return neighbours[x].size() != 0;
}

void process_input() {
    int parent, child;
    unsigned long parent_index, child_index;
    std::map<int, unsigned long>::iterator parent_iter, child_iter;
    while(file_stream >> parent >> child) {
        child_iter = renaming.find(child);
        if(child_iter == renaming.end()) {
            neighbours.push_back(std::vector<unsigned long>{});
            child_index = neighbours.size() - 1;
            renaming.insert(std::make_pair(child, child_index));
        }
        else child_index = child_iter->second;

        parent_iter = renaming.find(parent);
        if(parent_iter == renaming.end()) {
            neighbours.push_back(std::vector<unsigned long>{});
            parent_index = neighbours.size() - 1;
            renaming.insert(std::make_pair(parent, parent_index));
        }
        else parent_index = parent_iter->second;

        neighbours[parent_index].push_back(child_index);
    }
    n = neighbours.size();
    global_BC.resize(n, 0.);
}

long get_unprocessed() {
    {
        std::lock_guard<std::mutex> lock{next_mut};
        return (next_unprocessed >= n ? -1 : static_cast<long>(next_unprocessed++));
    }
}

bool has_unprocessed() {
    {
        std::lock_guard<std::mutex> lock{next_mut};
        return next_unprocessed < n;
    }
}

void process_graph() {
    if(!has_unprocessed()) return;
    std::vector<double> BC(n, 0.);
    std::vector<double> delta(n, 0.), sigma(n, 0.);
    std::vector<int> d(n, -1);
    std::vector<std::vector<unsigned long> > p{n, std::vector<unsigned long>{}};
    std::stack<unsigned long> reach_stack;
    std::queue<unsigned long> Q;
    while(true) {
        long s1 = get_unprocessed();
        if(s1 < 0) break;
        unsigned long s = static_cast<unsigned long>(s1);
 //       std::cout << s << " of " << n << std::endl;
        if(!has_neighbours(s)) continue;
        for(unsigned long i = 0; i < n; ++i) {
            delta[i] = sigma[i] = 0;
            d[i] = -1;
            p[i].clear();
        }

        sigma[s] = 1;
        d[s] = 0;
        Q.push(s);
        while(!Q.empty()) {
            unsigned long v = Q.front();
            Q.pop();
            reach_stack.push(v);
            for(unsigned long w = 0; w < neighbours[v].size(); ++w) {
                unsigned long neighbour = neighbours[v][w];
                if(d[neighbour] < 0) {
                    Q.push(neighbour);
                    d[neighbour] = d[v] + 1;
                }

                if(d[neighbour] == d[v] + 1) {
                    sigma[neighbour] += sigma[v];
                    p[neighbour].push_back(v);
                }
            }
        }

        while(!reach_stack.empty()) {
            unsigned long w = reach_stack.top();
            reach_stack.pop();
            for(auto v = p[w].begin(); v != p[w].end(); ++v)
                delta[*v] += ((sigma[*v] / sigma[w]) * (1. + delta[w]));

            if(w != s) BC[w] += delta[w];
        }

    }

    {
        std::lock_guard<std::mutex> lock{BC_mut};
        for (unsigned long i = 0; i < n; ++i) global_BC[i] += BC[i];
    }
}

void create_threads() {
    for(int i = 0; i < THREAD_AMOUNT; ++i) {
        thread_vector.push_back(std::thread(process_graph));
    }
    for(int i = 0; i < THREAD_AMOUNT; ++i) {
        thread_vector[i].join();
    }
}

int main(int argc, char* argv[]) {
    assert(argc == 4);
    THREAD_AMOUNT = atoi(argv[1]);
    INPUT_NAME = argv[2];
    OUTPUT_NAME = argv[3];
    file_stream.open(INPUT_NAME, std::ios::in);
    if(!file_stream.is_open()) return 1;

    process_input();

    create_threads();

    file_stream.close();
    file_stream.open(OUTPUT_NAME, std::ios::out);
    if(!file_stream.is_open()) return 1;

    for(auto& iter = renaming.begin(); iter != renaming.end(); ++iter) {
        if(has_neighbours(iter->second))
            file_stream << iter->first << " " << global_BC[iter->second] << std::endl;
    }

    file_stream.close();
    return 0;
}

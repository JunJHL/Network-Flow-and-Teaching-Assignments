#include "netflow.hpp"

#include <iostream>

int main()
{
    // The example below uses the flow network from slide #61 of
    // the graphs lecture slide deck, although the letters do not
    // correspond perfectly with 0 through 5, e.g. "a" from the slide
    // is not 0 here; it is 3. (This is to show you that the source
    // need not be the first vertex and the sink need not be the
    // last vertex.)
    std::vector<Edge> maxFlow = solveNetworkFlow({
            {3, 5, 5},
            {3, 0, 10},
            {0, 5, 5},
            {5, 4, 10},
            {0, 4, 3},
            {0, 1, 1},
            {4, 1, 20},
            {4, 2, 5},
            {1, 2, 7}
        },
        6);
    for (const Edge& edge : maxFlow)
        std::cout << edge.from << " -> " << edge.to
            << " (" << edge.weight << ")\n";
}

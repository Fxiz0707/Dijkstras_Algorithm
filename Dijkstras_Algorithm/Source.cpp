#include <iostream>
#include <vector>
#include <list> 

int g_Inf = 32000; 

using namespace std; 

struct Edge
{
	int mEndNode;
	int mWeight;

	Edge() : mEndNode(-1), mWeight(-1) {}; 
	Edge(int endNode, int weight) : mEndNode(endNode), mWeight(weight)
	{
	}
};

struct Graph
{
private:
	std::vector<std::list<Edge>> mData;

public: 
	const int mSize; 

	Graph(int size) : mSize(size), mData()
	{
		mData.resize(mSize); 
	}

	void add_edge(int u, int v, int weight)
	{
		mData[u].push_back(Edge(v, weight)); 
	}

	std::list<Edge> get_edges(int u)
	{
		return mData[u]; 
	}
};


/**
 * \brief Min-Heap, specifically to store unvisited nodes for Dijkstra's algorithm 
 */

struct MinHeap
{
private:
	std::vector<int> mData;
	std::vector<int> mPositions; // mPositions[i] holds the position in mData of node 'i'

public:

	/**
	 * \brief 
	 * \param size The size of the graph with which the min-heap is being used
	 */
	MinHeap(int size)
	{
		mPositions.resize(size); 
	}

	/**
	 * \brief 
	 * \param u index of the node that you are adding (the value for the heap)  
	 * \param distances[i] represents the current shortest distance to node 'i' (the key for the heap) 
	 */
	void insert(int u, int* distances)
	{
		mData.push_back(u);
		int i = mData.size() - 1;
		mPositions[u] = i; 

		// Bubbling UP 
		while (i != 0) 
		{
			if (distances[mData[i]] < distances[mData[(i - 1) / 2]]) // parent key > child key
			{
				// swap parent and child key 

				int temp = mData[(i - 1) / 2];
				mData[(i - 1) / 2] = mData[i];
				mData[i] = temp;

				mPositions[mData[(i - 1) / 2]] = i; 
				mPositions[mData[i]] = (i - 1)/2;

				i /= 2;
			}
			else break; 
		}
	}

	/**
	 * \brief To be called whenever edges are relaxed within Dijkstra's algorithm
	 * \param u index of the node which has been updated 
	 * \param distances distances[i] represents the current shortest distance to node 'i' (the key for the heap)
	 */
	void decrease_key(int u, int* distances) // It is assumed, within the distance array, only one key has been changed 
	{
		int i = mPositions[u];
		if (i == -1) return; 

		// Bubbling UP
		while (i != 0)
		{
			if (distances[mData[i]] < distances[mData[(i - 1) / 2]]) // parent key > child key
			{
				// swap parent and child key

				int temp = mData[(i - 1) / 2];
				mData[(i - 1) / 2] = mData[i];
				mData[i] = temp;

				mPositions[mData[(i - 1) / 2]] = i;
				mPositions[mData[i]] = (i - 1) / 2;

				i /= 2;
			}
			else break; 
		}
	}

	/**
	 * \brief Pop the un-visited node with the shortest distance from the source (the smallest key) 
	 * \param distances distances[i] represents the current shortest distance to node 'i' (the key for the heap)
	 * \return The id/index of the node popped 
	 */
	int extract_min(int* distances) 
	{
		int min = mData[0];
		mData[0] = mData[mData.size() - 1]; // swapping the first element in the array with the last, so that we can bubble down
		mData.pop_back();

		if (mData.size() == 0) return min; // There is nothing to bubble, the last node has been popped from the Heap 

		mPositions[mData[0]] = 0;
		mPositions[min] = -1; // This node is no longer in the Heap 

		int i = 0;
		// Bubbling DOWN 
		while (i < mData.size())
		{
			int smallerChildIndex = -1;
			if (distances[2 * i + 1] < distances[2 * i + 2]) // We need to swap with the smaller of the two child nodes 
			{
				smallerChildIndex = 2 * i + 1;
			}
			else smallerChildIndex = 2 * i + 2;

			if (distances[smallerChildIndex] < distances[i]) // parent key > smaller of the two child keys 
			{
				// swap parent and child key 

				int temp = mData[i];
				mData[i] = mData[smallerChildIndex];
				mData[smallerChildIndex] = temp;

				mPositions[mData[i]] = smallerChildIndex;
				mPositions[mData[smallerChildIndex]] = i;

				i = smallerChildIndex;
			}
			else break; 
		}

		return min; 
	}

	/**
	 * \brief Used as a stop condition to decide when to finish Dijkstra's Algorithm
	 * \return if the heap is empty
	 */
	bool is_empty()
	{
		return (mData.size() == 0); 
	}
};


/**
 * \brief Dijkstra's running in time complexity O((n + m) * log(n)) 
 * \param graph 
 */
void dijkstras_algorithm(Graph graph)
{
	int* shortest = new int[graph.mSize];
	shortest[0] = 0;

	for (int i = 1; i < graph.mSize; i++)
	{
		shortest[i] = g_Inf; 
	}

	MinHeap queue = MinHeap(graph.mSize);

	for (int i = 0; i < graph.mSize; i++)
	{
		queue.insert(i, shortest); 
	}

	while (!queue.is_empty()) // O((n + m) * log(n))
	{

		// extract_min is ran once for every node. Therefore adds time complexity O(n * log(n))
		int minIndex = queue.extract_min(shortest);

		std::list<Edge> edges = graph.get_edges(minIndex); 

		for (auto it = edges.begin(); it != edges.end(); it++)  
		{
			if (shortest[minIndex] + it->mWeight < shortest[it->mEndNode]) // relaxing the edge is successful 
			{
				shortest[it->mEndNode] = shortest[minIndex] + it->mWeight;

				// decrease_key is ran once for every single edge. Therefore adds time complexity O(m * log(n))
				queue.decrease_key(it->mEndNode, shortest); 
			}
		}
	}
}

 int main()
{
	Graph graph = Graph(6);
	graph.add_edge(0, 1, 5); 
	graph.add_edge(0, 2, 2); 
	graph.add_edge(0, 3, 6); 
	graph.add_edge(1, 4, 4); 
	graph.add_edge(2, 3, 2); 
	graph.add_edge(2, 5, 12); 
	graph.add_edge(3, 4, 4); 
	graph.add_edge(3, 5, 8); 
	graph.add_edge(4, 5, 3);

	graph.add_edge(1, 0, 5);
	graph.add_edge(2, 0, 2);
	graph.add_edge(3, 0, 6);
	graph.add_edge(4, 1, 4);
	graph.add_edge(3, 2, 2);
	graph.add_edge(5, 2, 12);
	graph.add_edge(4, 3, 4);
	graph.add_edge(5, 3, 8);
	graph.add_edge(5, 4, 3);

	dijkstras_algorithm(graph);
}


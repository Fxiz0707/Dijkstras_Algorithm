#include <iostream>
#include <vector>
#include <list> 

int g_Inf = 32000;

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
	std::vector<int> mTree;
	std::vector<int> mPositions; // mPositions[i] holds the position in mTree of node 'i' (decrease_key is the only function where mPositions is accessed)

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
		mTree.push_back(u);
		int i = mTree.size() - 1;
		mPositions[u] = i;

		// Bubbling UP 
		while (i != 0)
		{
			int parentKeyIndex = (i - 1) / 2;

			if (distances[mTree[i]] < distances[mTree[parentKeyIndex]]) // child key < parent key 
			{
				// swap parent and child key 
				int temp = mTree[parentKeyIndex];
				mTree[parentKeyIndex] = mTree[i];
				mTree[i] = temp;

				mPositions[mTree[parentKeyIndex]] = parentKeyIndex;
				mPositions[mTree[i]] = i;

				i = parentKeyIndex;
			}
			else
			{
				break;
			}
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

		while (i != 0)
		{
			int parentNodeIndex = (i - 1) / 2;

			if (distances[mTree[i]] < distances[mTree[parentNodeIndex]]) // child key < parent key 
			{
				// swap parent and child key
				int temp = mTree[parentNodeIndex];
				mTree[parentNodeIndex] = mTree[i];
				mTree[i] = temp;

				mPositions[mTree[parentNodeIndex]] = parentNodeIndex;
				mPositions[mTree[i]] = i;

				i = parentNodeIndex;
			}
			else
			{
				break;
			}
		}
	}

	/**
	 * \brief Pop the un-visited node with the shortest distance from the source (the smallest key)
	 * \param distances distances[i] represents the current shortest distance to node 'i' (the key for the heap)
	 * \return The id/index of the node popped
	 */
	int extract_min(int* distances)
	{
		int min = mTree[0];
		mTree[0] = mTree[mTree.size() - 1]; // swapping the first element in the array with the last, so that we can bubble down
		mTree.pop_back();

		if (mTree.size() == 0) return min; // There is nothing to bubble, the last node has been popped from the Heap 

		mPositions[mTree[0]] = 0;
		mPositions[min] = -1; // This node is no longer in the Heap

		int i = 0;

		// Bubbling DOWN 
		while (i < mTree.size())
		{
			int smallerChildIndex = -1;
			int leftChildIndex = 2 * i + 1;
			int rightChildIndex = 2 * i + 2;

			if (leftChildIndex >= mTree.size())
			{
				break;
			}

			if (rightChildIndex >= mTree.size())
			{
				smallerChildIndex = leftChildIndex;
			}
			else if (distances[mTree[leftChildIndex]] <= distances[mTree[rightChildIndex]])
			{
				smallerChildIndex = leftChildIndex;
			}
			else
			{
				smallerChildIndex = rightChildIndex;
			}

			if (distances[mTree[smallerChildIndex]] < distances[mTree[i]]) // child key < parent key 
			{
				// swap parent and child key 
				int temp = mTree[i];
				mTree[i] = mTree[smallerChildIndex];
				mTree[smallerChildIndex] = temp;

				mPositions[mTree[i]] = i;
				mPositions[mTree[smallerChildIndex]] = smallerChildIndex;

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
		return (mTree.size() == 0);
	}
};

/**
 * \brief Dijkstra's running in time complexity O((n + m) * log(n))
 * \param graph
 */
void dijkstras_algorithm(Graph graph)
{
	int* shortest = new int[graph.mSize]; // shorest[i] returns the current shortest distance from node 0 to node 'i'
	int* preds = new int[graph.mSize]; // preds[i] returns the predecessor to node 'i' in its current shortest path from node 0 
	preds[0] = 0;
	shortest[0] = 0;

	for (int i = 1; i < graph.mSize; i++)
	{
		shortest[i] = g_Inf;
		preds[i] = -1;
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
				preds[it->mEndNode] = minIndex;

				// decrease_key is ran once for every single edge. Therefore adds time complexity O(m * log(n))
				queue.decrease_key(it->mEndNode, shortest);
			}
		}
	}

	// Printing path from source to last node

	if (preds[graph.mSize - 1] == -1)
	{
		std::cout << "No path exists" << std::endl;
		return;
	}

	int node = graph.mSize - 1;
	while (node != 0)
	{
		std::cout << node << " <- ";
		node = preds[node];
	}
	std::cout << '0' << std::endl;
}

int main()
{
	Graph graph = Graph(7);

	graph.add_edge(0, 1, 6);
	graph.add_edge(0, 3, 10);
	graph.add_edge(0, 4, 18);
	graph.add_edge(1, 2, 12);
	graph.add_edge(1, 4, 9);
	graph.add_edge(3, 4, 14);
	graph.add_edge(3, 6, 22);
	graph.add_edge(4, 6, 16); 
	graph.add_edge(4, 2, 2); 
	graph.add_edge(2, 6, 13);
	graph.add_edge(2, 5, 27);
	graph.add_edge(6, 5, 8); 

	dijkstras_algorithm(graph);
}



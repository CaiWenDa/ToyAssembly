#pragma once
#include "overlap.h"
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/copy.hpp>
#include "boost/graph/properties.hpp"
#include "boost/format.hpp"


namespace cwd
{
	typedef struct AVertex {
		uint r;
		dnaPos_t SP;
		dnaPos_t EP;
		bool orient = true;

		friend ostream& operator<<(ostream& out, const AVertex& v)
		{
			out << v.r;
			return out;
		}
	} AVertex;

	typedef struct AEdge {
		enum Adj
		{
			HeadHead, HeadTail, TailHead, TailTail
		} adj;

		assemblyInfo_t ovl;
		double weight = 1.0;

		friend ostream& operator<<(ostream& out, const AEdge& e)
		{
			out << bitset<2>(e.adj);
			return out;
		}

	} AEdge;

	struct vertex_property_t {
		typedef boost::vertex_property_tag kind;
	};

	struct edge_property_t {
		typedef boost::edge_property_tag kind;
	};

	template <class WeightMap, class weight>
	class edge_writer
	{
	public:
		edge_writer(WeightMap w, weight n) : wm(w), wn(n) {}

		template <class Edge>
		void operator()(ostream& out, const Edge& e) const {
			out << "[adj=" << bitset<2>(wm[e]) << ", weight=" << wn[e] << "]";
		}
	private:
		WeightMap wm;
		weight wn;
	};

	template <class WeightMap, class weight>
	inline edge_writer<WeightMap, weight> make_edge_writer(WeightMap w, weight n)
	{
		return edge_writer<WeightMap, weight>(w, n);
	}

	using AGraph = boost::adjacency_list<boost::mapS, boost::vecS, boost::bidirectionalS, boost::property<vertex_property_t, AVertex>, AEdge>;
	using SubGraph = boost::adjacency_list<boost::mapS, boost::vecS, boost::undirectedS, boost::property<vertex_property_t, AVertex>, AEdge>;
	using ComponentGraph = boost::filtered_graph<AGraph, function<bool(AGraph::edge_descriptor)>, function<bool(AGraph::vertex_descriptor)> >;

	using vertex_iterator = boost::graph_traits<AGraph>::vertex_iterator;
	using edge_iterator = boost::graph_traits<AGraph>::edge_iterator;
	using edge_descriptor = boost::graph_traits<AGraph>::edge_descriptor;
	using vertex_descriptor = boost::graph_traits<AGraph>::vertex_descriptor;
	using out_edge_iterator = boost::graph_traits<AGraph>::out_edge_iterator;
	using in_edge_iterator = boost::graph_traits<AGraph>::in_edge_iterator;


	void generateContig(const cwd::AGraph& g, std::vector<vertex_descriptor>& p, const cwd::seqData_t& seq, std::ofstream& seqOut, int id, int& totalLen);
	uint heuristic(edge_descriptor & p, const cwd::AGraph& g, boost::adj_list_edge_property_map<boost::bidirectional_tag, cwd::AEdge::Adj, const cwd::AEdge::Adj&, cwd::vertex_descriptor, const cwd::AEdge, cwd::AEdge::Adj cwd::AEdge::*>& adj_prop, boost::adj_list_edge_property_map<boost::bidirectional_tag, cwd::assemblyInfo_t, const cwd::assemblyInfo_t&, cwd::vertex_descriptor, const cwd::AEdge, cwd::assemblyInfo_t cwd::AEdge::*>& ovl_prop);
	vector<vector<vertex_descriptor>> findPath(const AGraph& g, const seqData_t& seq, ofstream& outAssembly);
	seqan::Dna5String concatReads(const seqan::Dna5String& pre, uint r2, const seqan::Dna5String& R2, const assemblyInfo_t& ovl, cwd::AEdge::Adj adj, bool shouldRev, bool toPrev);
	seqan::Dna5String concatReadsDirect(const seqan::Dna5String& pre, const seqan::Dna5String& R2);
	void createOverlapGraph(seqData_t& seq, size_t block1, size_t block2);
	void assembler(const seqData_t& seq, ofstream& seqOut);
	void readOverlapGraph(const string& graphName, vector<AGraph>& v_g);
	void connected_components_subgraphs(AGraph const& g);
}
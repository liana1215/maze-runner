#ifndef MAZE_APP_HPP
#define MAZE_APP_HPP

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <SDL/SDL.h>

namespace MazeApp
{
    using namespace std;

    class MazeBuilder {
    friend class MazeSolver;
    friend class MazeDrawer;

    public:
        MazeBuilder(uint32_t);
        ~MazeBuilder()=default;

        struct Cell {
            uint32_t row;
            uint32_t col;
        };

        struct Channel {
            Cell src;
            Cell dst;
        };

        struct  Vertex {
            int32_t number; 
            Cell cell;      
        };

        void buildMaze();
        void runAldousBroder();

    private:
        Vertex getVertex(uint32_t, vector<vector<Vertex>>);
        Cell getNeighborLocation(uint32_t, uint32_t, uint32_t);
        bool checkNeighborExists(uint32_t, uint32_t, uint32_t);

    private:
        vector<vector<vector<uint32_t>>> maze;
        vector<Channel> channels;
    };


    class MazeSolver {
    friend class MazeDrawer;

    public:
        MazeSolver(uint32_t);
        ~MazeSolver()=default;

        struct Path {
            stack<MazeBuilder::Vertex> p;
        };

        Path solveMazeBFS(MazeBuilder &);

    private:
        void getGrid();
        
    private:
        uint32_t m_dim;
        vector<vector<uint32_t>> grid;
    };


    class MazeDrawer {
    public:
        MazeDrawer(uint32_t, uint32_t);
        ~MazeDrawer()=default;

        uint32_t drawMaze(MazeBuilder &, MazeSolver&, vector<MazeBuilder::Vertex>);

    private:
        uint32_t checkLoad(SDL_Surface* );
        vector<uint32_t> getCellType(MazeBuilder&, MazeSolver&, uint32_t);

    private:
        uint32_t bmp_size;
        uint32_t n_cell_template;
    };
}  

#endif

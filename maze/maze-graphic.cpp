#include "MazeApp.hpp"

using namespace MazeApp;

int 
main(int argc, char** argv)
{
    if (argc == 1) {
        std::cerr << "Input size of desired maze as parameter.\n";
        std::exit(1);
    }
    else {
        auto maze_dim = stoi(argv[1]);
        
        MazeBuilder builder = MazeBuilder(maze_dim);
        MazeSolver solver = MazeSolver(maze_dim);
        MazeDrawer drawer = MazeDrawer(32, 16);

        builder.buildMaze();
        builder.runAldousBroder();
        
        MazeSolver::Path solved = solver.solveMazeBFS(builder);
        std::vector<MazeBuilder::Vertex> shortest_path;
        auto stack_length = solved.p.size();    

        std::cout << std::endl;
        std::cout << "Maze solved: ";
        std::cout << "GOAL: ";
        for (auto i = 0; i < stack_length; i++) { 
            MazeBuilder::Vertex v = solved.p.top(); 
            shortest_path.push_back(v);
            std::cout << v.number << " <-";
            solved.p.pop();
        }
        std::cout << "START" << std::endl;

        drawer.drawMaze(builder, solver, shortest_path);
        return 0;
    }
}


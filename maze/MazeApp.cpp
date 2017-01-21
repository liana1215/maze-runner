#include "MazeApp.hpp"

#define PRINTBMP  // dont define if you want sample outputs of maze graphics
#define GIFBMP  // dont define if you want sample outputs for gif of maze graphics 

using namespace MazeApp;

MazeBuilder::MazeBuilder(uint32_t dim = 10)
    : maze(dim, vector<vector<uint32_t>>(dim)) {}


void
MazeBuilder::buildMaze()
{
    auto m_dim = maze.size();
    uint32_t direction[] = {0, 1, 2, 3};   // {up, right, down, left}
    
    for (auto r = 0; r < m_dim; r++) 
        for (auto c = 0; c < m_dim; c++) 
            for (auto d : direction) {
                if (r == 0 && c == 0) {
                    maze[r][c].push_back((d == 1 || d == 2));
                }
                else if (r == m_dim - 1 && c == m_dim - 1) {
                    maze[r][c].push_back((d == 0 || d == 3));
                }
                else if (r == 0 && c == m_dim - 1) {
                    maze[r][c].push_back((d == 2 || d == 3));
                }
                else if (r == m_dim - 1 && c == 0) {
                    maze[r][c].push_back((d == 0 || d == 1));  
                }
                else if (r == 0) {
                    maze[r][c].push_back((d == 1 || d == 2 || d == 3));   
                }
                else if (c == 0) {
                    maze[r][c].push_back((d == 0 || d == 1 || d == 2)); 
                }
                else if (r == m_dim - 1) {
                    maze[r][c].push_back((d == 0 || d == 1 || d == 3)); 
                }
                else if (c == m_dim - 1) {
                    maze[r][c].push_back((d == 0 || d == 2 || d == 3)); 
                }
                else {
                    maze[r][c].push_back(1);
                }
            }
}
    

void        
MazeBuilder::runAldousBroder()
{
    random_device rd;
    mt19937 gen(rd());
    uint32_t maze_size = maze.size();
    uniform_int_distribution<> dist(0,maze_size-1); 
    uniform_int_distribution<> neighbor(0,3);
    
    vector<vector<uint32_t>> visited(maze_size, vector<uint32_t>(maze_size));
        
    uint32_t r = 0;
    uint32_t c = 0;
    visited[r][c] = 1;
    uint32_t remaining = maze_size * maze_size - 1;
    uint32_t count = 0;
    
    while (remaining > 0) {
        uint32_t n_index = neighbor(gen);
        Cell neighbor_cell;

        while(!checkNeighborExists(r, c, n_index)) {
            n_index = neighbor(gen);
        }
        neighbor_cell = getNeighborLocation(r, c, n_index);

        if (!visited[neighbor_cell.row][neighbor_cell.col]) {
            Channel channel;
            channel.src.row = r;
            channel.src.col = c;
            channel.dst.row = neighbor_cell.row;
            channel.dst.col = neighbor_cell.col;
            channels.push_back(channel);       //add channel to vector of channels
            visited[neighbor_cell.row][neighbor_cell.col] = 1;
            remaining--;
        }
        r = neighbor_cell.row;
        c = neighbor_cell.col;

        count++;
    }
    
    for (auto i = 0; i < channels.size(); i++) {
        cout << channels[i].src.row << "," << channels[i].src.col << " -> "
             << channels[i].dst.row << "," << channels[i].dst.col << endl;
    }
    cout << "Steps taken to construct maze: " << count << endl;
}


MazeBuilder::Cell
MazeBuilder::getNeighborLocation(uint32_t r, uint32_t c, uint32_t index)
{
    auto r_n = r;
    auto c_n = c;
    Cell n;

    switch (index) {
        case 0: r_n = r - 1; break;
        case 1: c_n = c + 1; break;
        case 2: r_n = r + 1; break;
        case 3: c_n = c - 1; break;
    }
    n.row = r_n;
    n.col = c_n;
    return n;
}


bool
MazeBuilder::checkNeighborExists(uint32_t r, uint32_t c, uint32_t n_index)
{
    return maze[r][c][n_index];
}


MazeBuilder::Vertex
MazeBuilder::getVertex(uint32_t number, vector<vector<MazeBuilder::Vertex>> graph)
{
    for (auto v : graph)
        for (auto d : v) 
            if (number == d.number) {
                return d;
            }
    throw invalid_argument("start vertex not valid"); // vertex not found
}


MazeSolver::MazeSolver(uint32_t dim = 10) 
    : grid(dim, vector<uint32_t>(dim)), m_dim(dim) {}


void
MazeSolver::getGrid()
{
    uint32_t n = 0;
    for (int r = 0; r < m_dim; r++) {
        for (int c = 0; c < m_dim; c++) {
            grid[r][c] = n;
            n++;
        }
    }
}


MazeSolver::Path
MazeSolver::solveMazeBFS(MazeBuilder &builder) 
{
    Path path_to_goal; // return value
    queue<Path> paths; // queue container

    uint32_t size = m_dim * m_dim;
    uint32_t start_vertex = 0;      // starting cell
    uint32_t goal = size-1;         // last cell

    vector<string> c(size);         // color
    vector<uint32_t> d(size);       // distance
    vector<uint32_t> test(size);    // test for existance
    vector<int32_t> p(size);        // predecessor

    getGrid();

    vector<vector<MazeBuilder::Vertex>> graph(size);
    for (auto c: builder.channels) {
        MazeBuilder::Vertex one, two;
        one.cell = c.src;
        one.number = grid[c.src.row][c.src.col];
        two.cell = c.dst;
        two.number = grid[c.dst.row][c.dst.col];

        test[one.number] = 1;
        test[two.number] = 1;
        graph[one.number].push_back(two);
        graph[two.number].push_back(one);
    }

    for (int i = 0; i < graph.size(); i++) 
        if (i != start_vertex) {
            if (test[i]) {
                c[i] = "WHITE";
                d[i] = 10000;
                p[i] = -1;
            } 
            else {
                c[i] = "BLACK";
            }
        } 
        else {
            c[i] = "GRAY";
            d[i] = 0;
            p[i] = -1;
        }

    Path start;     //initialize first path
    MazeBuilder::Vertex s;

    try {
        s = builder.getVertex(start_vertex, graph);
    } 
    catch (const exception& e) {
        cout << "ERROR: " << e.what() << endl;
    }

    start.p.push(s);            //push vertex number to stack 
    paths.push(start);          //enque path

    while (!paths.empty()) {
        Path u = paths.front();
        MazeBuilder::Vertex check = u.p.top();

        if (check.number == goal) {
            path_to_goal = u;
        } 
        else {
            for (int i = 0; i < graph[check.number].size(); i++) {
                MazeBuilder::Vertex v = graph[check.number][i];
                if (c[v.number] == "WHITE") {
                    Path path_copy = u;
                    path_copy.p.push(v);
                    paths.push(path_copy);
                    c[v.number] = "GRAY";
                    d[v.number] = d[check.number] + 1;
                    p[v.number] = check.number;
                }
            }
        }
        paths.pop();
        c[check.number] = "BLACK";
    }
    return path_to_goal;
}


MazeDrawer::MazeDrawer(uint32_t sz, uint32_t cell)
    : bmp_size(sz), n_cell_template(cell) {}


uint32_t
MazeDrawer::drawMaze(MazeBuilder& builder, MazeSolver& solver,
                     vector<MazeBuilder::Vertex> shortest_path)
{
    vector<MazeBuilder::Channel> channels = builder.channels;
    uint32_t maze_dim = solver.m_dim;
    
    SDL_Surface *screen;
    SDL_Surface *monster;
    SDL_Surface *cell_template;
    SDL_Surface *start;
    SDL_Surface *goal;
    SDL_Surface *trail;

    SDL_Rect src, dest;
    string monster_bmp = "./static/m32.bmp";
    string cell_bmp = "./static/c32.bmp";
    string start_bmp = "./static/s32.bmp";
    string goal_bmp = "./static/g32.bmp";
    string trail_bmp = "./static/trail.bmp";

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cout << "Unable to initialize SDL: " <<  SDL_GetError()
                  << endl;
        return 1;
    }
    atexit(SDL_Quit);
    uint32_t window_dim = maze_dim;
    if (maze_dim < 5) { window_dim = 5; }

    screen = SDL_SetVideoMode(window_dim*bmp_size,
                              window_dim*bmp_size,
                              bmp_size,
                              0);

    if (screen == NULL) {
        cout << "Unable to set video mode: " << SDL_GetError()
                  << endl;
        return 1;
    }
    monster = SDL_LoadBMP(monster_bmp.c_str()); checkLoad(monster);
    cell_template = SDL_LoadBMP(cell_bmp.c_str()); checkLoad(cell_template);
    start = SDL_LoadBMP(start_bmp.c_str()); checkLoad(start);
    goal  = SDL_LoadBMP(goal_bmp.c_str()); checkLoad(goal);
    trail = SDL_LoadBMP(trail_bmp.c_str()); checkLoad(trail);

    vector<uint32_t>types = getCellType(builder, solver, maze_dim * maze_dim);

    for (int r = 0; r < maze_dim; r++) 
        for (int c = 0; c < maze_dim; c++) {
            uint32_t type = types[solver.grid[r][c]];
            src.x = type * bmp_size;
            src.y = 0;
            src.w = (cell_template->w) / n_cell_template;
            src.h = cell_template->h;
            dest.x = c * bmp_size;
            dest.y = r * bmp_size;
            dest.w = (cell_template->w) / n_cell_template;
            dest.h = cell_template->h;
            SDL_BlitSurface(cell_template, &src, screen, &dest);
            SDL_UpdateRect(screen, 0, 0, 0, 0);
        }
    
    // Place monster at start
    uint32_t colorkey = SDL_MapRGB(monster->format, 0, 0, 0);
    SDL_SetColorKey(monster, SDL_SRCCOLORKEY, colorkey);
    src.x = 0;
    src.y = 0;
    src.w = monster->w;
    src.h = monster->h;
    dest.x = 0;
    dest.y = 0;
    dest.w = monster->w;
    dest.h = monster->h;
    SDL_BlitSurface(monster, &src, screen, &dest);
    SDL_UpdateRect(screen, 0, 0, 0, 0);

    // Place goal icon
    SDL_SetColorKey(goal, SDL_SRCCOLORKEY, colorkey);
    src.x = 0;
    src.y = 0;
    src.w = goal->w;
    src.h = goal->h;
    dest.x = (maze_dim - 1) * bmp_size;
    dest.y = (maze_dim - 1) * bmp_size;
    dest.w = monster->w;
    dest.h = monster->h;
    SDL_BlitSurface(goal, &src, screen, &dest);
    SDL_UpdateRect(screen, 0, 0, 0, 0);

    // Output solution path
    SDL_SetColorKey(trail, SDL_SRCCOLORKEY, colorkey);
    for (int i = shortest_path.size()-2; i > 0; i--) {
        // start from 2nd cell, as 1st cell has monster icon.
        MazeBuilder::Vertex v = shortest_path[i];
        uint32_t r = v.cell.row;
        uint32_t c = v.cell.col;
        src.x = 0;
        src.y = 0;
        src.w = (trail->w);
        src.h = trail->h;
        dest.x = c*bmp_size;
        dest.y = r*bmp_size;
        dest.w = trail->w;
        dest.h = trail->h;
        SDL_BlitSurface(trail, &src, screen, &dest);
        SDL_UpdateRect(screen, 0, 0, 0, 0);
        SDL_Delay(500);

        // save surface to bmp as samples
#ifndef PRINTBMP
        if (i == (shortest_path.size()-2)) {
            SDL_SaveBMP(screen, "start");
        }
        if (i == (shortest_path.size()/2)) {
            SDL_SaveBMP(screen, "middle");
        }
        if (i == 1) {
            SDL_SaveBMP(screen, "end");
        }
#endif

#ifndef GIFBMP
        string gif = ("gif_" + to_string((shortest_path.size()-2)-i) + ".bmp");
        SDL_SaveBMP(screen, gif.c_str());
#endif
    }

    SDL_Delay(5000);

    /* Free the memory that was allocated to the bitmap. */
    SDL_FreeSurface(monster);
    SDL_FreeSurface(cell_template);
    SDL_FreeSurface(start);
    SDL_FreeSurface(goal);
    return 0;
}


uint32_t
MazeDrawer::checkLoad(SDL_Surface* image)
{
    if (image == NULL) {
        cout << "Unable to load bitmap.\n"; 
        return 1;
    }
}


vector<uint32_t>
MazeDrawer::getCellType(MazeBuilder& builder, MazeSolver& solver, uint32_t n_cells) 
{
    vector<vector<uint32_t>> neighbors(n_cells, vector<uint32_t>(4));
    // where neighbors [0,1,2,3] (u,r,d,l)
    vector<uint32_t> types(n_cells);
    
    for (auto c : builder.channels) {
        uint32_t cell_s = solver.grid[c.src.row][c.src.col];
        uint32_t cell_d = solver.grid[c.dst.row][c.dst.col];

        if ((c.dst.row == c.src.row) && (c.dst.col == (c.src.col + 1))) {
            neighbors[cell_s][1] = 1;
            neighbors[cell_d][3] = 1;
        } 
        else if ((c.dst.row == c.src.row) && (c.dst.col == (c.src.col - 1))) {  
            neighbors[cell_s][3] = 1;
            neighbors[cell_d][1] = 1;
        } 
        else if ((c.dst.col == c.src.col) && (c.dst.row == (c.src.row + 1))) {
            neighbors[cell_s][2] = 1;
            neighbors[cell_d][0] = 1;
        } 
        else if ((c.dst.col == c.src.col) && (c.dst.row == (c.src.row - 1))) {
            neighbors[cell_s][0] = 1;
            neighbors[cell_d][2] = 1;
        }
    }

    for (auto &n : neighbors) 
    {
        auto i = &n - &neighbors[0];
        if      (n[0]==0 && n[1]==1 && n[2]==1 && n[3]==1) {types[i] = 0;}
        else if (n[0]==1 && n[1]==1 && n[2]==0 && n[3]==1) {types[i] = 1;}
        else if (n[0]==1 && n[1]==0 && n[2]==1 && n[3]==1) {types[i] = 2;}
        else if (n[0]==1 && n[1]==1 && n[2]==1 && n[3]==0) {types[i] = 3;}
        else if (n[0]==0 && n[1]==1 && n[2]==1 && n[3]==0) {types[i] = 4;} 
        else if (n[0]==0 && n[1]==0 && n[2]==1 && n[3]==1) {types[i] = 5;} 
        else if (n[0]==1 && n[1]==1 && n[2]==0 && n[3]==0) {types[i] = 6;} 
        else if (n[0]==1 && n[1]==0 && n[2]==0 && n[3]==1) {types[i] = 7;} 
        else if (n[0]==0 && n[1]==1 && n[2]==0 && n[3]==1) {types[i] = 8;} 
        else if (n[0]==1 && n[1]==0 && n[2]==1 && n[3]==0) {types[i] = 9;} 
        else if (n[0]==1 && n[1]==1 && n[2]==1 && n[3]==1) {types[i] = 10;} 
        else if (n[0]==1 && n[1]==0 && n[2]==0 && n[3]==0) {types[i] = 11;} 
        else if (n[0]==0 && n[1]==1 && n[2]==0 && n[3]==0) {types[i] = 12;} 
        else if (n[0]==0 && n[1]==0 && n[2]==1 && n[3]==0) {types[i] = 13;} 
        else if (n[0]==0 && n[1]==0 && n[2]==0 && n[3]==1) {types[i] = 14;} 
        else                                               {types[i] = 15;}
    }
    return types;
}

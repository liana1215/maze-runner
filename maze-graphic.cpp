#include <SDL/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>

//#define PRINTBMP  // define if you want sample outputs of maze graphics
//#define GIFBMP  // define if you want sample outputs for gif of maze graphics

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

struct Path {
    std::stack<Vertex> p;
};

 
void                    buildMaze(std::vector<std::vector<std::vector<uint32_t>>> &);
Cell                    getNeigborLocation(uint32_t, uint32_t, uint32_t);
bool                    neighborExists(std::vector<std::vector<std::vector<uint32_t>>>,
                                       uint32_t, uint32_t, uint32_t);

std::vector<Channel>    runAldousBroder(std::vector<std::vector<std::vector<uint32_t>>> &);
void                    getGrid(std::vector<std::vector<uint32_t>> &, uint32_t); 
Vertex                  getVertex(uint32_t, std::vector<std::vector<Vertex>>);
Path                    solveMazeBFS(std::vector<Channel>, uint32_t); 

uint32_t                checkLoad(SDL_Surface *);
std::vector<uint32_t>   getCellType(std::vector<std::vector<uint32_t>>, std::vector<Channel>, uint32_t);
uint32_t                drawMaze(std::vector<Channel>, std::vector<Vertex>, uint32_t,uint32_t, uint32_t); 


int main()
{
    uint32_t maze_dim = 10;
    std::cout << "Input size of maze: ";
    std::cin >> maze_dim;
    std::vector<std::vector<std::vector<uint32_t>>> maze(maze_dim, 
                                                         std::vector<std::vector<uint32_t>>(maze_dim));
    buildMaze(maze);
    std::vector<Channel> channels = runAldousBroder(maze);
    Path solved = solveMazeBFS(channels, maze_dim);
    std::vector<Vertex> shortest_path;

    uint32_t stack_length = solved.p.size();    
    std::cout << std::endl;
    std::cout << "Maze solved: ";
    std::cout << "GOAL: ";
    for (int i = 0; i < stack_length; i++) { 
        Vertex v = solved.p.top(); 
        shortest_path.push_back(v);
        std::cout << v.number << " <-";
        solved.p.pop();
    }
    std::cout << "START" << std::endl;
    
    drawMaze(channels, shortest_path, maze_dim, 32, 16);
    return 0;
}


void 
buildMaze(std::vector<std::vector<std::vector<uint32_t>>> &maze) 
{
    uint32_t maze_dim    = maze.size();
    uint32_t direction[] = {0,1,2,3}; // {Up, Right, Down, Left}

    for (uint32_t r = 0; r < maze_dim; r++) {
        for (uint32_t c = 0; c < maze_dim; c++) {
            for (auto d : direction) {
                if (r==0 && c==0) {
                    maze[r][c].push_back((d==1 || d==2)); 
                } 
                else if (r==maze_dim-1 && c==maze_dim-1) {
                    maze[r][c].push_back((d==0 || d==3));
                }
                else if (r==0 && c==maze_dim-1) {
                    maze[r][c].push_back((d==2 || d==3));
                } 
                else if (r==maze_dim-1 && c==0) {
                    maze[r][c].push_back((d==0 || d==1));
                } 
                else if (r==0) {
                    maze[r][c].push_back((d==1 || d==2 || d==3)); 
                }
                else if (c==0) {
                    maze[r][c].push_back((d==0 || d==1 || d==2));
                }
                else if (r==maze_dim-1) {
                    maze[r][c].push_back((d==0 || d==1 || d==3)); 
                }
                else if (c==maze_dim-1) {
                    maze[r][c].push_back((d==0 || d==2 || d==3));
                }
                else {
                    maze[r][c].push_back(1);
                }
             } 
        }
    }
}
    

Cell 
getNeigborLocation(uint32_t r, uint32_t c, uint32_t index)
{
    uint32_t r_n = r;
    uint32_t c_n = c;
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
checkNeighborExists(std::vector<std::vector<std::vector<uint32_t>>> const maze, 
                    uint32_t r, uint32_t c, uint32_t n_index) 
{
    return maze[r][c][n_index];
}


std::vector<Channel>
runAldousBroder(std::vector<std::vector<std::vector<uint32_t>>> & maze) 
{
    std::random_device rd;
    std::mt19937 gen(rd());
    uint32_t maze_size = maze.size();
    std::uniform_int_distribution<> dist(0,maze_size-1);
    std::uniform_int_distribution<> neighbor(0,3);
    std::vector<std::vector<uint32_t>> visited(maze_size,
                                               std::vector<uint32_t>(maze_size));
    std::vector<Channel> channels;

    uint32_t r         = 0;     // dist(gen);
    uint32_t c         = 0;     // dist(gen);
    visited[r][c]      = 1;     // mark initial point as visited     
    uint32_t remaining = maze_size * maze_size - 1;
    int count          = 0;

    while (remaining > 0) {
        uint32_t n_index = neighbor(gen);
        Cell neighbor_cell;
        while(!checkNeighborExists(maze, r, c, n_index)) {
            n_index = neighbor(gen);
        }
        neighbor_cell = getNeigborLocation(r,c,n_index);
        if (!visited[neighbor_cell.row][neighbor_cell.col]) {
            Channel channel;
            channel.src.row = r;
            channel.src.col = c;
            channel.dst.row = neighbor_cell.row;
            channel.dst.col = neighbor_cell.col;
            channels.push_back(channel);                 //add channel to vector of channels
           
            visited[neighbor_cell.row][neighbor_cell.col] = 1;
            remaining--;
        }
        r = neighbor_cell.row;                      
        c = neighbor_cell.col;

            count++;
    }
    for (int i = 0; i < channels.size(); i++) {
        std::cout << channels[i].src.row << "," << channels[i].src.col << " -> "
                  << channels[i].dst.row << "," << channels[i].dst.col << std::endl;
    }
    std::cout << "Steps taken to construct maze: " << count << std::endl;
    return channels;
} 


Vertex
getVertex(uint32_t number, std::vector<std::vector<Vertex>> graph)
{
    for (auto v: graph) {
        for (auto d: v) {
            if (number == d.number) {
                return d;
            }
        }
    }
    throw std::invalid_argument("start vertex not valid");        // vertex not found
}


void
getGrid(std::vector<std::vector<uint32_t>> &grid, uint32_t dim) 
{
    uint32_t n = 0;            // counter variable
    for (int r = 0; r < dim; r++) {
        for (int c = 0; c < dim; c++) {
            grid[r][c] = n; 
            n++;
        }
    }
}
    

Path
solveMazeBFS(std::vector<Channel> channels, uint32_t dim) 
{
    Path                     path_to_goal; // return value
    std::queue<Path>         paths;        // queue container

    uint32_t size          = dim * dim;
    uint32_t start_vertex  = 0;            // starting cell
    uint32_t goal          = size-1;       // last cell

    std::vector<std::string> c(size);      // color
    std::vector<uint32_t>    d(size);      // distance
    std::vector<uint32_t>    test(size);   // test for existance
    std::vector<int32_t>     p(size);      // predecessor

    std::vector<std::vector<uint32_t>> grid(dim, std::vector<uint32_t>(dim));
    getGrid(grid, dim);

    std::vector<std::vector<Vertex>> graph(size);
    for (auto c : channels) {
        Vertex one, two;
        one.cell   = c.src;
        one.number = grid[c.src.row][c.src.col];
        two.cell   = c.dst;
        two.number = grid[c.dst.row][c.dst.col];

        test[one.number] = 1;
        test[two.number] = 1;
        graph[one.number].push_back(two);
        graph[two.number].push_back(one);
    }

    for (int i = 0; i < graph.size(); i++) {
        if (i != start_vertex) {
            if (test[i]) {
                c[i] = "WHITE";
                d[i] = 10000;
                p[i] = -1;
            } else {
                c[i] = "BLACK";
            }
        } else {
            c[i] = "GRAY";
            d[i] = 0;
            p[i] = -1;
        }    
    }

    Path start;                 //initialize first path
    Vertex s;

    try {
        s = getVertex(start_vertex, graph);
    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
    }

    start.p.push(s);            //push vertex number to stack 
    paths.push(start);          //enque path

    while (!paths.empty()) {
        Path u       = paths.front();
        Vertex check = u.p.top();

        if (check.number == goal) {
            path_to_goal = u;
        } else {
            for (int i = 0; i < graph[check.number].size(); i++) {
                Vertex v = graph[check.number][i];
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


uint32_t 
checkLoad(SDL_Surface *image) 
{
    if (image == NULL) {
        std::cout << "Unable to load bitmap." 
                  << std::endl;
        return 1;
    }
}


std::vector<uint32_t>
getCellType(std::vector<std::vector<uint32_t>> grid, 
            std::vector<Channel> channels, 
            uint32_t n_cells) {
    std::vector<std::vector<uint32_t>> neighbors(n_cells,
                                                 std::vector<uint32_t>(4));
    // where neighbors [0,1,2,3] (u,r,d,l)
    std::vector<uint32_t> types(n_cells);
    
    for (auto c: channels) {
        uint32_t cell_s = grid[c.src.row][c.src.col];
        uint32_t cell_d = grid[c.dst.row][c.dst.col];

        if ((c.dst.row == c.src.row) && 
            (c.dst.col == (c.src.col + 1))) {
            neighbors[cell_s][1] = 1;
            neighbors[cell_d][3] = 1;
        } else if ((c.dst.row == c.src.row) && 
                   (c.dst.col == (c.src.col - 1))) {  
            neighbors[cell_s][3] = 1;
            neighbors[cell_d][1] = 1;
        } else if ((c.dst.col == c.src.col) && 
                   (c.dst.row == (c.src.row + 1))) {
            neighbors[cell_s][2] = 1;
            neighbors[cell_d][0] = 1;
        } else if ((c.dst.col == c.src.col) && 
                   (c.dst.row == (c.src.row - 1))) {
            neighbors[cell_s][0] = 1;
            neighbors[cell_d][2] = 1;
        }
    }

    for (auto &n: neighbors) {
        uint32_t i = &n - &neighbors[0];
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


uint32_t
drawMaze(std::vector<Channel> channels, std::vector<Vertex> shortest_path,
         uint32_t maze_dim, uint32_t bmp_size, uint32_t n_cell_template)
{
    SDL_Surface *screen;
    SDL_Surface *monster;
    SDL_Surface *cell_template;
    SDL_Surface *start;
    SDL_Surface *goal;
    SDL_Surface *trail;

    SDL_Rect src, dest;
    std::string monster_bmp = "m32.bmp";
    std::string cell_bmp    = "c32.bmp";
    std::string start_bmp   = "s32.bmp";
    std::string goal_bmp    = "g32.bmp";
    std::string trail_bmp   = "trail.bmp";

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cout << "Unable to initialize SDL: " <<  SDL_GetError()
                  << std::endl;
        return 1;
    }
    atexit(SDL_Quit);
    uint32_t window_dim = maze_dim;
    if (maze_dim < 5) {window_dim = 5;}

    screen = SDL_SetVideoMode(window_dim*bmp_size,
                              window_dim*bmp_size,
                              bmp_size,
                              0);

    if (screen == NULL) {
        std::cout << "Unable to set video mode: " << SDL_GetError()
                  << std::endl;
        return 1;
    }
    monster       = SDL_LoadBMP(monster_bmp.c_str()); checkLoad(monster);
    cell_template = SDL_LoadBMP(cell_bmp.c_str()); checkLoad(cell_template);
    start         = SDL_LoadBMP(start_bmp.c_str()); checkLoad(start);
    goal          = SDL_LoadBMP(goal_bmp.c_str()); checkLoad(goal);
    trail         = SDL_LoadBMP(trail_bmp.c_str()); checkLoad(trail);
    std::vector<std::vector<uint32_t>> grid(maze_dim,
                                            std::vector<uint32_t>(maze_dim));
    getGrid(grid, maze_dim);

    std::vector<uint32_t>types = getCellType(grid, channels, maze_dim*maze_dim);

    for (int r = 0; r < maze_dim; r++) {
        for (int c = 0; c < maze_dim; c++) {
            uint32_t type = types[grid[r][c]];
            src.x   = type*bmp_size;
            src.y   = 0;
            src.w   = (cell_template->w)/n_cell_template;
            src.h   = cell_template->h;
            dest.x  = c*bmp_size;
            dest.y  = r*bmp_size;
            dest.w  = (cell_template->w)/n_cell_template;
            dest.h  = cell_template->h;
            SDL_BlitSurface(cell_template, &src, screen, &dest);
            SDL_UpdateRect(screen, 0, 0, 0, 0);
        }
    }
    // Place monster at start
    uint32_t colorkey = SDL_MapRGB(monster->format, 0, 0, 0);
    SDL_SetColorKey(monster, SDL_SRCCOLORKEY, colorkey);
    src.x   = 0;
    src.y   = 0;
    src.w   = monster->w;
    src.h   = monster->h;
    dest.x  = 0;
    dest.y  = 0;
    dest.w  = monster->w;
    dest.h  = monster->h;
    SDL_BlitSurface(monster, &src, screen, &dest);
    SDL_UpdateRect(screen, 0, 0, 0, 0);

    // Place goal icon
    SDL_SetColorKey(goal, SDL_SRCCOLORKEY, colorkey);
    src.x   = 0;
    src.y   = 0;
    src.w   = goal->w;
    src.h   = goal->h;
    dest.x  = (maze_dim-1)*bmp_size;
    dest.y  = (maze_dim-1)*bmp_size;
    dest.w  = monster->w;
    dest.h  = monster->h;
    SDL_BlitSurface(goal, &src, screen, &dest);
    SDL_UpdateRect(screen, 0, 0, 0, 0);

    // Output solution path
    SDL_SetColorKey(trail, SDL_SRCCOLORKEY, colorkey);
    for (int i = shortest_path.size()-2; i > 0; i--) {
        // start from 2nd cell, as 1st cell has monster icon.
        Vertex v = shortest_path[i];
        uint32_t r = v.cell.row;
        uint32_t c = v.cell.col;
        src.x   = 0;
        src.y   = 0;
        src.w   = (trail->w);
        src.h   = trail->h;
        dest.x  = c*bmp_size;
        dest.y  = r*bmp_size;
        dest.w  = trail->w;
        dest.h  = trail->h;
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
        std::string gif = ("gif_" +
                            std::to_string((shortest_path.size()-2)-i)
                            + ".bmp");
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

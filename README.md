# Maze-Runner program

Automatically generates a m x m perfect maze, where the perfect maze is defined
as a maze with one solution, path to goal.

---

Usage:
```
$ make
$ ./mazeapp 10
```
The program takes in one parameter, m, which is used to define the dimension of
the maze. 

The program builds and solves the maze and outputs the path to goal, steps taken
to solve the maze, and a graphic that dynamically changes as the "monster"
approaches the goal.

Currently, the solver uses the aldous broder algorithm, which is clearly
inefficient. Implementation of more efficient algorithms to follow.

Sample output:
```
Steps taken to construct maze: 1022

Maze solved: GOAL: 99 <-98 <-97 <-96 <-86 <-85 <-95 <-94 <-84 <-74 <-64 <-54
<-44 <-43 <-33 <-34 <-24 <-14 <-15 <-5 <-4 <-3 <-2 <-1 <-0 <-START
```

![](https://github.com/liana1215/maze-runner/blob/master/maze/static/www.GIFCreator.me_RPu34o.gif)

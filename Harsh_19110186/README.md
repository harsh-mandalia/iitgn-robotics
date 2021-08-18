General instructions to run the code:
1. Each task has one or more than one python file that you need to run to see the
animation.
2. Press “q” while “window” is selected to close the window. Or press ctrl+c while the
terminal is selected to force quit.
3. Feel free to change the parameters in the \_\_init__(self) function.
4. Task1.py involved kinematic approach where calculated q1 and q2 are directly animated,
but task1_dynamics.py involves dynamics approach where tau1 and tau2 (torques) are
calculated. <br />
    > In task1.py the bot will follow the mouse pointer if it's under "window". <br />
    > In task1_dynamics.py first draw the path by draging the mouse on the "window" then the bot will follow the given path.
5. In task2.py, self.wall1 and self.wall2 are points of two ends of the wall. The wall is being
shown in animation between those two points only but actually, the wall is taken as an
infinite line for the 2R manipulator.
6. In task3.py self.x0 and self.y0 are the mean position of the spring(bot).
7. In task4.py in line 47 (mybot.trace(5)) change the value 5 to something else (integer) to
change the speed at with the simulation draws the workspace.

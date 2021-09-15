path = "d:\Program Files (x86)\Graphviz2.38\bin\";%path%

del ex6_task1.gif
del ex6_task1b.gif
del ex6_task2.gif
del ex6_task2_new.gif
dot -Tgif ex6_task1.dot -oex6_task1.gif
dot -Tgif ex6_task1b.dot -oex6_task1b.gif
dot -Tgif ex6_task2.dot -oex6_task2.gif
dot -Tgif ex6_task2_new.dot -oex6_task2_new.gif

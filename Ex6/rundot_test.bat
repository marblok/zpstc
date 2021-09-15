path = "d:\Program Files (x86)\Graphviz2.38\bin\";%path%

del ex6_task2_test.gif
dot -Tgif ex6_task2_test.dot -oex6_task2_test.gif
del ex6_task2.gif
dot -Tgif ex6_task2.dot -oex6_task2.gif
del Seria2Parallel-test.gif
dot -Tgif Seria2Parallel-test.dot -oSeria2Parallel-test.gif
del SymbolMapper-test.gif
dot -Tgif SymbolMapper-test.dot -oSymbolMapper-test.gif


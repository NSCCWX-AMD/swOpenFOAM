find . -name "*.c" |xargs cat|grep -v ^$|wc -l
find . -name "*.h" |xargs cat|grep -v ^$|wc -l
find . -name "*.hpp" |xargs cat|grep -v ^$|wc -l
find . -name "*.cpp" |xargs cat|grep -v ^$|wc -l
find . -name "*.C" |xargs cat|grep -v ^$|wc -l
find . -name "*.H" |xargs cat|grep -v ^$|wc -l
find . -name "*.sh" |xargs cat|grep -v ^$|wc -l

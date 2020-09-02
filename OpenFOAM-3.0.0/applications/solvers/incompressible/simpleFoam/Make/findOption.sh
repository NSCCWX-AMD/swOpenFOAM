#!/bin/sh
#************** get options libs *********************
find_key="files"
option_key="options"
special_key="libOpenFOAM"

cd ./Make
cur_path=$PWD
find_path=${WM_PROJECT_DIR}/src
cat options | grep -w "EXE_LIBS_NEW" 2>&1 > /dev/null
if [ $? -eq 0 ]; then
    rm -f options.copy
else
    cp options options.bak
    cat options | grep -w "EXE_LIBS" 2>&1 > /dev/null
    if [ $? -eq 0 ]; then
        cat options | grep '\-l' > options.copy 
        sed -e 's/EXE_LIBS =//g' -e "s/ *-l/lib/g" -e "s/ //g" -e 's/\\//g' options.copy > options.tmp
        strtmp=$(cat options.tmp); 
        strtmp=$(echo ${strtmp} | sed -e 's/ /:/g' )
        echo ${strtmp}
        rm -f options.copy options.tmp
        
        #*************** trunc strtmp to a array *****************
        index=0
        OLDIFS=$IFS
        IFS=":"
        for item in ${strtmp}; do
        a[${index}]="$item"
        a_append[${index}]="$item"
        index=`expr ${index} + 1`
        done
        IFS=$OLDIFS
        
        function find_func(){
        
            local lib_tmp=$1
            local b
            local i
            
            cd ${find_path}
            echo $PWD
            echo $find_key
            #cd $find_tmp
            find ./ -name "${find_key}" > findfile.tmp
            files=$(cat findfile.tmp)
            for Macrofile in ${files}; do
			    echo ${Macrofile}
                grep -w "${lib_tmp}" ${Macrofile} 2>&1 > /dev/null
                if [ $? -eq 0 ]; then
        #            echo ${Macrofile} 
                     filepath=$(echo ${Macrofile}|sed "s/$(basename ${Macrofile})//g")
        #            echo ${filepath}
                     optionfile=${filepath}"options"
                     echo ${optionfile}
            #        echo "11111111111111111"
                     cat ${optionfile} | grep " \-l" 2>&1 > /dev/null
                     if [ $? -eq 0 ]; then
                         str_b=$(cat ${optionfile} | grep '\ -l' | sed -e 's/LIB_LIBS =//g' -e "s/ *-l/lib/g" -e "s/ //g" -e 's/\\//g')
					     str_b=$(echo ${str_b} |  sed -e 's/ /:/g') 
			echo str_b: $str_b
                        echo ${lib_tmp} depends  ${str_b} 
            #            echo "22222222222222222"
                #3333333333333333333333333333333
                        jndex=0
                        OLDIFS=$IFS
                        IFS=":"
                        for item in ${str_b}; do
                            b[$jndex]="$item"
                            jndex=`expr ${jndex} + 1`
                        done
                        IFS=$OLDIFS
                #3333333333333333333333333333333
                    
                
            #    	     echo "33333333333333333"
                #4444444444444444444444444444444
                         for ((i=0;i<${#b[@]};i++)); do
                             echo $strtmp | grep -w "${b[$i]}" 2>&1 > /dev/null
                             #echo $?
                             if [ $? -ne 0 ]; then
                                 if [ ${b[$i]} = ${special_key} ]; then
				     echo Do nothing for libOpenFOAM
                                     #strtmp=$strtmp":"${special_key}":libPstream"
                                     #find_func "libPstream"
                                 else
            #                        echo "4444444444444444444",$strtmp               
                                     strtmp=$strtmp":"${b[$i]}
            #                        echo "555555555555555555",$strtmp
                                     find_func ${b[$i]}  
                                 fi
                             fi
                         done
            #4444444444444444444444444444444  
                    fi
                fi  
            
            #2222222222222222222222222222222222222222222222222222222
            done
            rm -f findfile.tmp
        }
        
        for ((i=0;i<${#a[@]};i++)); do
            echo "a[" $i "]:" ${a[$i]}
            cd ${cur_path}
            if [ ${a[$i]} = ${special_key} ]; then
                echo Do nothing for libOpenFOAM
                #strtmp=${strtmp}":libPstream"
                #find_func "libPstream"
            else
                echo "find_func" ${a[$i]}
                find_func ${a[$i]}
            fi
        done
        strtmp="EXE_LIBS_NEW = "$(echo ${strtmp} | sed -e "s/:/ /g" -e "s/lib/-l/g")
        #echo ${strtmp} 
        cd ${cur_path};echo ${strtmp} > options.new
        #sed -i '$a ########## add by hexiang ############' options
        sed -i '$ r options.new' options
        rm -f options.new
    fi
fi
cat options

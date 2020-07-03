## This script generates a seris of directires in which there is a relax.in
## Author: Yuewen FANG; Contact: fyuewen@gmail.com; Date: 2020 July 2nd
####original cell parameter#####
a=[   2.512269491   1.450459467   4.397688755
  -2.512269491   1.450459467   4.397688755
   0.000000000  -2.900918933   4.397688755]

for i in 93:107
    ratio = i*0.01
    round_ratio = round(ratio, digits=2) # keep two digits
    rm(string(round_ratio), recursive=true)
    mkdir(string(round_ratio))
    filename = string(string(round_ratio),"/relax.in")
    println(filename)
    fout = open(filename, "w")
    ###the part before the CELLPARAMETER###
    f = open("relax.in") 
    lines = readlines(f)
    for i in lines
        if startswith(i, "   2.512269")
            break
        else
                write(fout, i)
                write(fout, '\n')
        end
    end
    
    new_cell_parameter = cbrt(ratio)*a
    for row in eachrow(new_cell_parameter)
        for arr_element in row
            write(fout, string(arr_element))
            write(fout, ' ')
#             print(arr_element)
#             print(' ')
        end
        write(fout, '\n')
    end

    ###the part after ATOMIC_POSITIONS###
    f = open("relax.in") 
    lines = readlines(f)
    line_n = 0
    for i in lines
#        global line_n
        line_n = line_n + 1
        if startswith(i, "ATOMIC_POSITIONS")
            break
        else
#                 write(fout, i)
#                 write(fout, '\n')
        end
    end  
#    print(line_n)
            
    
    line_nn = 0
    for i in lines
#        local line_nn
        line_nn = line_nn+1
        if line_nn >= line_n
            write(fout, i)
            write(fout, '\n')
        end
    end
    
    close(fout)
end


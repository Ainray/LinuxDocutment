########################################################################
#                                                                      #
# author:        Ainray                                                #
#                                                                      #
# email:         wwzhang0421@163.com                                   #
#                                                                      #
# file:          bash                                                  #
#                                                                      #
# created:       2017-12-24 20:20:50                                   #
#                                                                      #
# modified:      2017-12-24 23:22:30                                   #
#                                                                      #
# introduction:  Someting about Bash. This a not a working script,     #
#the first line just is for using vim syntax colors.                   #
#                                                                      #
# license:       Public license.                                       #
#                                                                      #
########################################################################
help::
    help set

debug::
    bash -x script > log

    debug parts in a script by,
        set -x
            debug part
        set +x

    change debug mode
        set -f  ( or set -o noglob) disable file name generation 
        set -v  ( or set -o verbose) prints shell input lines
        set -x  ( or set -o xtrace) print command traces before executing command

    use echo to print message

Quote::
        Quotes preserve literal meaning of special character, that is, it stops
    Bash from intepretating it and maybe expanding it.
        Example,
            ls [Vv]*    #[,],* are special characters
            ls '[Vv]*'  #now nothing in quotes is special, ls try to find file called '[Vv]*' exactly.
        Some certain programs and utilites reintepret or expand special characters in a quoted string.
    This because, Bash first preserve and parse then into "the program" but removing the surrouding 
    quotes, so the protected characters again exposed to "the program" which, like Bash, interpret and 
    expand these special characters.
        soft/double quotes,
            + $,`,\(follewed by $ and `) in double quotes is still special
            + double quotes stop word spliting
            Example,
                var=""
                COMMAND $var $var $var          # execute COMMAND with no arguments
                COMMAND "$var" "$var" "$var"    # execute COMMAND with 3 empty arguments
                COMMAND "$var $var $var"        # execute COMMAND with 1 arguments (two spaces)
        strong/single quotes,
            + $ is literal
            + \ is literal but except in cases of $'\033'or $'\x27' 

jobs::
        kill some hanged job,
            job & { sleep ${TIMEOUT}; eval 'kill -9 $!' &> /dev/null;} #this maybe some bug
            count=0
            job & { 
                while ((count < TIMEOUT));do
                    eval '[ ! -d "/proc/$lastjob" ]  && ((count=TIMEOUT))'
                    last=$!  #record its pid
                    ((count++))
                    sleep 1
                done
            eval ' [ -d "/proc/$lastjob" ] && kill -15 $lastjob' 
            }

read::
        Shell builtin.
        refer keydetect function

regular expression::

test::
        (()) and let test arithmetic, if artithmetic is zero, the return is 1; otherwise, it is 0
            Example,
                var=-2 && ((var+=2 )) && echo $var  # because the second arthimetic expression gives 0, so its
                                                # exit status is non-zero, Bash thinks some error occurs and
                                                # do not print $var
                let "num = (( 200|11))"         # $? gives 0, number gives 203
        if test any command,
            if cmd 
            then
            else
            fi
            Example
                if cmp a b &> /dev/null
                then echo "Files a and b are identical"
                else echo "Files a and b differ"
                fi
                word=linux
                if echo "${word}" | grep -q "inu"
                then echo "inu found in $word"
                else echo "inu not found in $word"
                fi
        test; /usr/bin/test; []; and /usr/bin/[; are the equivalence
            if test -z "$1";
            if /usr/bin/test -z "$1";
            if [ -z "$1" ];
            if /usr/bin[ -z "$1";
        [[]], arithmetic evaluation of octal/hexadecimal constants takes place automatically within a [[]]
            decimal=15
            octol=017
            hex=0x0f
            if [ "$decimal" -eq "$octal" ]; 
            if [[ "$decimal" -eq "$octal" ]]; 

source::
        check a script is sourced or not, by 
        [[ $_ != $0 ]] && ISSOURCE=1 || ISSOURCE=0

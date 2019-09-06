special characters::
    
    space, \s

number of occurrences::
    \+, one or more
    *, zero or more

upper::
lower::
    sed 's/[A-Z]/\l&/g'
    sed 's/[A-Z]/\u&/g'

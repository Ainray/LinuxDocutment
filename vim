function::
   1. definition syntax
      :fu[nction][!] {name}([arguments]) [range] [abort] [dict] [closure]
      function, vim key word
      !, force redefine funtion with the same name as a existed function
      {name}, first char must be caplitalized
      [arguments], conveyed as a pointer-like way
Tricks::
  command-line window:
    q:(in normal mode) to open it
  convert matched word into upper-case ones:
    %s/iptables/\U&/g
  copy current command:
    :let hits=GetMatches(654,723,'shift',1,1)  
    ":p
  delete matched lines
    :g/^ *$/d
    :g!/^ *$/d  or :v/^ *$/d
    :v/error\|varn\|fail/d
  execute current yanked command:
    :@"
  one-line multiple replacement
    call feedkeys(repeat("yn",5))|s/\$/\\>\$/gc
  scrolling synchronously:
    set scb (in each window)
    set scb! (toggle on/off)
  vertical to horizontal or vice versa:
    Ctrl-w t Ctrl-w H  (h->v)
      t, make first top left, H 
      H, make current to full-height at far left
    Ctrl-w t Ctrl-w K (v-h)
      K, make current to full-width at the very top
     

Spellcheck::
  set spelllang=en
  set spell "to start spell
  ]s "jump to next misspelled word
  z= "see suggestions
  zg "add misspelled word under cursor to your personal spellfile, maybe ~/.vim/spell/en.utf-8.add
  zG "ignore the misspelled word under cursor for this ssesion
  zw "to mark correctly spelled word under cursor as a misspeeling

regular::

    *   any number
    \+  one or more
    \w  [a-zA-z]

    example:
        /^ *local\+      "matching leading local

Syntastic::
    disable Syntastic auto check by 
        :SyntasticToggleMode
    do no show flake8 quick fix
        let g:flake8_show_quickfix=0
    then recheck by press <F7> or call 
        call flake8#Flake8()

word separator::
    temporarily allow - as a in word for * check word by,
        set iskeyword=+-

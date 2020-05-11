move::
    1. key word hyperlink: Ctrl+] and Ctrl+T
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
  search current word or character
    :yw
    /<C-R><C-R>"
    :yl
  search case sensitive
    /\Cpi
  search ignore case
    /\cpi
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
   change width and height of current window
      :res +5 
      :res -5
      :vertical resize +5

Spellcheck::
  set spelllang=en
  set spell "to start spell
  ]s "jump to next misspelled word
  z= "see suggestions
  zg "add misspelled word under cursor to your personal spellfile, maybe ~/.vim/spell/en.utf-8.add
  zG "ignore the misspelled word under cursor for this ssesion
  zw "to mark correctly spelled word under cursor as a misspeeling
regular::
zero::
count::
    *   any number
    \+  one or more
    \w  [a-zA-z]
    example:
        /^ *local\+      "matching leading local
Special characters::
    code of current character,
        ga  " ascii
        g8  " utf8
    enter by code,
        <C+v> 64    "@
        <C+v> u64   "d
        <C+k> <c1><c2>
        <C+k> 12  " ½
        @
    digital graphs,
        :digraphs  " check codes
        <C+k> 12   " enter by digital graphs
        
Syntastic::
    do no show flake8 quick fix let g:flake8_show_quickfix=0
    then recheck by press <F7> or call 
        call flake8#Flake8()
word separator::
    temporarily allow - as a in word for * check word by,
        set iskeyword=+-
    yw
    yl
     
help::
helptags::
tags::
    helptags ~/.vim/bundle/vim-latex/doc/
display::
list::
    display line end character "$" and so on
        set list
count::
    characters,
        :%s/./&//gn  
    words,
        :%s/\i\+/&/gn
    lines,
        :%s/^//n
    pattern,
        :%s/pattern/&/gn  
    pattern as a word,
        :%s/\<the\>/&/gn
quantifier::
zero::
more::
one::
n::
    zero or more, as few as possible,
        \{-}
    one or more
        \{-n,m}
    at least n or more
        \{-n,}
    one or at most m
        \{-,m}

plugin::
taglist:: source code brower
    + dependence on ctags, install it by 
        sudo apt-get install exuberant-ctags
    + dependence on file type in .vimrc
        filetype on
    + vimrc general tags configuration
        set tags+=~/.vim/tags/cpp
        set tags+=~/.vim/tags/gl
        set tags+=~/.vim/tags/linuxheader
        set tags+=./.tags
    + vimrc configration
        let Tlist_Auto_Open = 1
        let Tlist_Auto_Update = 1
        let Tlist_Close_On_Select = 1
        let Tlist_Compact_Format = 1
        let Tlist_Ctags_Cmd = '/usr/local/bin/ctags'  "only necessary if it is not included in PATH
        let Tlist_Display_Prototype = 1
        let Tlist_Exist_OnlyWindow = 1
        let Tlist_File_Fold_Auto_Close = 1
    + commands
        :TlistOpen
        :TlistToggle

    + current tag, type P on current tag
    
    + toggle, in normal 
        map <silent> <F10> :TlistToggle<CR>


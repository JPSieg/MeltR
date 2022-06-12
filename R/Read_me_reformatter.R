list.files()

Con = file("README.md")
Lines = readLines(Con)
close(Con)

#Remove duplicates for dark mode

Remove.index = c()

for (i in 1:length(Lines)){
  if (length(grep("#gh-dark-mode-only", Lines[i])) > 0){
    Remove.index = c(Remove.index, i)
  }
}

Lines = Lines[-Remove.index]

for (i in 1:length(Lines)){

  if (length(grep("img src", Lines[i])) > 0){
    Lines[i] = gsub("\"", "", Lines[i])
    Lines[i] = sub(as.character("<img src= https://render.githubusercontent.com/render/math?math={"), "$$", Lines[i], fixed = TRUE)
    Lines[i] = gsub("}#gh-light-mode-only>", "$$", Lines[i], fixed = TRUE)
    Lines[i] = gsub("%2b", "+", Lines[i], fixed = TRUE)
  }else{
    Lines[i] = gsub("<sup>", "\textsuperscript{", Lines[i], fixed = TRUE)
    Lines[i] = gsub("</sup>", "}", Lines[i], fixed = TRUE)
    Lines[i] = gsub("<sub>", "\textsubscript{", Lines[i], fixed = TRUE)
    Lines[i] = gsub("</sub>", "}", Lines[i], fixed = TRUE)
    Lines[i] = as.character(Lines[i])
  }
}

Con = file("Manual.md")
writeLines(Lines, Con)
close(Con)


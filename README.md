# PW-Parallel-Brandes-Algorithm
Concurrent programming assignment #2

**Task description in Polish**

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<!-- saved from url=(0082)http://www.mimuw.edu.pl/~janowska/PW-LAB/Zadania/Zadanie2/zadanie-brandes-cpp.html -->
<html xmlns="http://www.w3.org/1999/xhtml" lang="pl" xml:lang="pl"><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<title>Zadanie zaliczeniowe CPP - algorytm Brandesa</title>
<!-- 2016-12-26 Mon 21:06 -->

<meta name="generator" content="Org-mode">
<meta name="author" content="Krzysztof Rządca">
<style type="text/css">
 <!--/*--><![CDATA[/*><!--*/
  .title  { text-align: center; }
  .todo   { font-family: monospace; color: red; }
  .done   { color: green; }
  .tag    { background-color: #eee; font-family: monospace;
            padding: 2px; font-size: 80%; font-weight: normal; }
  .timestamp { color: #bebebe; }
  .timestamp-kwd { color: #5f9ea0; }
  .right  { margin-left: auto; margin-right: 0px;  text-align: right; }
  .left   { margin-left: 0px;  margin-right: auto; text-align: left; }
  .center { margin-left: auto; margin-right: auto; text-align: center; }
  .underline { text-decoration: underline; }
  #postamble p, #preamble p { font-size: 90%; margin: .2em; }
  p.verse { margin-left: 3%; }
  pre {
    border: 1px solid #ccc;
    box-shadow: 3px 3px 3px #eee;
    padding: 8pt;
    font-family: monospace;
    overflow: auto;
    margin: 1.2em;
  }
  pre.src {
    position: relative;
    overflow: visible;
    padding-top: 1.2em;
  }
  pre.src:before {
    display: none;
    position: absolute;
    background-color: white;
    top: -10px;
    right: 10px;
    padding: 3px;
    border: 1px solid black;
  }
  pre.src:hover:before { display: inline;}
  pre.src-sh:before    { content: 'sh'; }
  pre.src-bash:before  { content: 'sh'; }
  pre.src-emacs-lisp:before { content: 'Emacs Lisp'; }
  pre.src-R:before     { content: 'R'; }
  pre.src-perl:before  { content: 'Perl'; }
  pre.src-java:before  { content: 'Java'; }
  pre.src-sql:before   { content: 'SQL'; }

  table { border-collapse:collapse; }
  caption.t-above { caption-side: top; }
  caption.t-bottom { caption-side: bottom; }
  td, th { vertical-align:top;  }
  th.right  { text-align: center;  }
  th.left   { text-align: center;   }
  th.center { text-align: center; }
  td.right  { text-align: right;  }
  td.left   { text-align: left;   }
  td.center { text-align: center; }
  dt { font-weight: bold; }
  .footpara:nth-child(2) { display: inline; }
  .footpara { display: block; }
  .footdef  { margin-bottom: 1em; }
  .figure { padding: 1em; }
  .figure p { text-align: center; }
  .inlinetask {
    padding: 10px;
    border: 2px solid gray;
    margin: 10px;
    background: #ffffcc;
  }
  #org-div-home-and-up
   { text-align: right; font-size: 70%; white-space: nowrap; }
  textarea { overflow-x: auto; }
  .linenr { font-size: smaller }
  .code-highlighted { background-color: #ffff00; }
  .org-info-js_info-navigation { border-style: none; }
  #org-info-js_console-label
    { font-size: 10px; font-weight: bold; white-space: nowrap; }
  .org-info-js_search-highlight
    { background-color: #ffff00; color: #000000; font-weight: bold; }
  /*]]>*/-->
</style>
<link href="./Zadanie zaliczeniowe CPP - algorytm Brandesa_files/cpp-theme.css" rel="stylesheet">
<script type="text/javascript">
/*
@licstart  The following is the entire license notice for the
JavaScript code in this tag.

Copyright (C) 2012-2013 Free Software Foundation, Inc.

The JavaScript code in this tag is free software: you can
redistribute it and/or modify it under the terms of the GNU
General Public License (GNU GPL) as published by the Free Software
Foundation, either version 3 of the License, or (at your option)
any later version.  The code is distributed WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU GPL for more details.

As additional permission under GNU GPL version 3 section 7, you
may distribute non-source (e.g., minimized or compacted) forms of
that code without the copy of the GNU GPL normally required by
section 4, provided you include this license notice and a URL
through which recipients can access the Corresponding Source.


@licend  The above is the entire license notice
for the JavaScript code in this tag.
*/
<!--/*--><![CDATA[/*><!--*/
 function CodeHighlightOn(elem, id)
 {
   var target = document.getElementById(id);
   if(null != target) {
     elem.cacheClassElem = elem.className;
     elem.cacheClassTarget = target.className;
     target.className = "code-highlighted";
     elem.className   = "code-highlighted";
   }
 }
 function CodeHighlightOff(elem, id)
 {
   var target = document.getElementById(id);
   if(elem.cacheClassElem)
     elem.className = elem.cacheClassElem;
   if(elem.cacheClassTarget)
     target.className = elem.cacheClassTarget;
 }
/*]]>*///-->
</script>
<script type="text/javascript" src="./Zadanie zaliczeniowe CPP - algorytm Brandesa_files/MathJax.js"></script>
<script type="text/javascript">
<!--/*--><![CDATA[/*><!--*/
    MathJax.Hub.Config({
        // Only one of the two following lines, depending on user settings
        // First allows browser-native MathML display, second forces HTML/CSS
        //  config: ["MMLorHTML.js"], jax: ["input/TeX"],
            jax: ["input/TeX", "output/HTML-CSS"],
        extensions: ["tex2jax.js","TeX/AMSmath.js","TeX/AMSsymbols.js",
                     "TeX/noUndefined.js"],
        tex2jax: {
            inlineMath: [ ["\\(","\\)"] ],
            displayMath: [ ['$$','$$'], ["\\[","\\]"], ["\\begin{displaymath}","\\end{displaymath}"] ],
            skipTags: ["script","noscript","style","textarea","pre","code"],
            ignoreClass: "tex2jax_ignore",
            processEscapes: false,
            processEnvironments: true,
            preview: "TeX"
        },
        showProcessingMessages: true,
        displayAlign: "center",
        displayIndent: "2em",

        "HTML-CSS": {
             scale: 100,
             availableFonts: ["STIX","TeX"],
             preferredFont: "TeX",
             webFont: "TeX",
             imageFont: "TeX",
             showMathMenu: true,
        },
        MMLorHTML: {
             prefer: {
                 MSIE:    "MML",
                 Firefox: "MML",
                 Opera:   "HTML",
                 other:   "HTML"
             }
        }
    });
/*]]>*///-->
</script>
<style type="text/css">#MathJax_Message {position: fixed; left: 1em; bottom: 1.5em; background-color: #E6E6E6; border: 1px solid #959595; margin: 0px; padding: 2px 8px; z-index: 102; color: black; font-size: 80%; width: auto; white-space: nowrap}
#MathJax_MSIE_Frame {position: absolute; top: 0; left: 0; width: 0px; z-index: 101; border: 0px; margin: 0px; padding: 0px}
.MathJax_Error {color: #CC0000; font-style: italic}
</style><style type="text/css">.MathJax_Display {text-align: center; margin: 1em 0em; position: relative; display: block; width: 100%}
.MathJax .merror {background-color: #FFFF88; color: #CC0000; border: 1px solid #CC0000; padding: 1px 3px; font-family: serif; font-style: normal; font-size: 90%}
.MathJax_Preview {color: #888888}
#MathJax_Tooltip {background-color: InfoBackground; color: InfoText; border: 1px solid black; box-shadow: 2px 2px 5px #AAAAAA; -webkit-box-shadow: 2px 2px 5px #AAAAAA; -moz-box-shadow: 2px 2px 5px #AAAAAA; filter: progid:DXImageTransform.Microsoft.dropshadow(OffX=2, OffY=2, Color='gray', Positive='true'); padding: 3px 4px; position: absolute; left: 0; top: 0; width: auto; height: auto; display: none}
.MathJax {display: inline; font-family: serif; font-style: normal; font-weight: normal; line-height: normal; font-size: 100%; font-size-adjust: none; text-indent: 0; text-align: left; text-transform: none; letter-spacing: normal; word-spacing: normal; word-wrap: normal; white-space: nowrap; float: none; border: 0; padding: 0; margin: 0}
.MathJax img, .MathJax nobr, .MathJax a {border: 0; padding: 0; margin: 0; max-width: none; max-height: none; vertical-align: 0; line-height: normal; text-decoration: none}
img.MathJax_strut {border: 0 !important; padding: 0 !important; margin: 0 !important; vertical-align: 0 !important}
.MathJax span {display: inline; border: 0; padding: 0; margin: 0; vertical-align: 0; line-height: normal; text-decoration: none}
.MathJax .MathJax_HitBox {cursor: text}
#MathJax_Tooltip * {filter: none; background: transparent}
#MathJax_getScales {font-family: STIXGeneral,STIXNonUnicode,STIXSizeOneSym}
.MathJax .math span {font-family: STIXGeneral,STIXNonUnicode,STIXSizeOneSym}
</style></head>
<body><div style="visibility: hidden; overflow: hidden; height: 1px; position: absolute; top: 0px;"><div id="MathJax_Hidden"></div></div><div id="MathJax_Message" style="display: none;"></div>
<div id="content">
<h1 class="title">Zadanie zaliczeniowe CPP - algorytm Brandesa</h1>
<div id="table-of-contents">
<h2>Spis treści</h2>
<div id="text-table-of-contents">
<ul>
<li><a href="http://www.mimuw.edu.pl/~janowska/PW-LAB/Zadania/Zadanie2/zadanie-brandes-cpp.html#sec-1">1. Wprowadzenie</a></li>
<li><a href="http://www.mimuw.edu.pl/~janowska/PW-LAB/Zadania/Zadanie2/zadanie-brandes-cpp.html#sec-2">2. Algorytm Brandesa</a></li>
<li><a href="http://www.mimuw.edu.pl/~janowska/PW-LAB/Zadania/Zadanie2/zadanie-brandes-cpp.html#sec-3">3. Wejście i wyjście</a></li>
<li><a href="http://www.mimuw.edu.pl/~janowska/PW-LAB/Zadania/Zadanie2/zadanie-brandes-cpp.html#sec-4">4. Forma oddania zadania</a></li>
<li><a href="http://www.mimuw.edu.pl/~janowska/PW-LAB/Zadania/Zadanie2/zadanie-brandes-cpp.html#sec-5">5. Kryteria oceny</a></li>
<li><a href="http://www.mimuw.edu.pl/~janowska/PW-LAB/Zadania/Zadanie2/zadanie-brandes-cpp.html#sec-6">6. Dodatkowe materiały</a></li>
<li><a href="http://www.mimuw.edu.pl/~janowska/PW-LAB/Zadania/Zadanie2/zadanie-brandes-cpp.html#sec-7">7. FAQ</a></li>
</ul>
</div>
</div>
<p>
Termin zadania: 9 stycznia 2017
</p>

<p>
Modyfikacje treści: 
</p>
<ul class="org-ul">
<li>2016/12/26 FAQ
</li>
<li>2016/12/23 FAQ
</li>
<li>2016/12/13 FAQ
</li>
<li>2016/12/09 ogłoszenie zadania
</li>
</ul>


<div id="outline-container-sec-1" class="outline-2">
<h2 id="sec-1"><span class="section-number-2">1</span> Wprowadzenie</h2>
<div class="outline-text-2" id="text-1">
<p>
Analiza sieci polega m.in. na identyfikacji <span class="underline">ważnych</span>
wierzchołków. Naiwnym sposobem wyrażania ważności wierzchołka jest
np. zliczenie sąsiadów—im więcej sąsiadów tym ważniejszy jest
wierzchołek dla sieci. Liczenie sąsiadów nie identyfikuje jednak
ważnych z punktu widzenia spójności sieci wierzchołków-pośredników
(np. wierzchołka łączącego dwie grupy niepołączonych ze sobą
wierzchołków). Zaproponowano wiele alternatywnych miar,
m.in. closeness centrality czy eigenvector centrality; tematem zadania
jest obliczenie pośrednictwa (betweenness).
</p>

<p>
Pośrednictwo węzła <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-1"><span style="display: inline-block; position: relative; width: 0.472em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.915em 1000em 2.727em -0.51em);"><span class="mrow" id="MathJax-Span-2"><span class="mi" id="MathJax-Span-3" style="font-family: STIXGeneral; font-style: italic;">v</span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.487em; vertical-align: -0.019em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-1">v</script> określa jaka część najkrótszych ścieżek w sieci
przechodzi przez ten węzeł: <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-4"><span style="display: inline-block; position: relative; width: 8.614em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.534em; left: 0em; clip: rect(1.333em 1000em 3.168em -0.539em);"><span class="mrow" id="MathJax-Span-5"><span class="mi" id="MathJax-Span-6" style="font-family: STIXGeneral; font-style: italic;">B</span><span class="mi" id="MathJax-Span-7" style="font-family: STIXGeneral; font-style: italic;">C</span><span class="mo" id="MathJax-Span-8" style="font-family: STIXGeneral;">(</span><span class="mi" id="MathJax-Span-9" style="font-family: STIXGeneral; font-style: italic;">v</span><span class="mo" id="MathJax-Span-10" style="font-family: STIXGeneral;">)</span><span class="mo" id="MathJax-Span-11" style="font-family: STIXGeneral; padding-left: 0.313em;">=</span><span class="texatom" id="MathJax-Span-12" style="padding-left: 0.313em;"><span class="mrow" id="MathJax-Span-13"><span class="munderover" id="MathJax-Span-14"><span style="display: inline-block; position: relative; width: 2.805em; height: 0px;"><span style="position: absolute; top: -2.534em; left: 0em; clip: rect(1.596em 1000em 2.972em -0.473em);"><span class="mo" id="MathJax-Span-15" style="font-family: STIXGeneral; vertical-align: -0.002em;">∑</span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; top: -2.238em; left: 0.96em;"><span class="texatom" id="MathJax-Span-16"><span class="mrow" id="MathJax-Span-17"><span class="mi" id="MathJax-Span-18" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">s</span><span class="mo" id="MathJax-Span-19" style="font-size: 70.7%; font-family: STIXGeneral;">≠</span><span class="mi" id="MathJax-Span-20" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">v</span><span class="mo" id="MathJax-Span-21" style="font-size: 70.7%; font-family: STIXGeneral;">≠</span><span class="mi" id="MathJax-Span-22" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">t</span></span></span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span></span></span><span class="mfrac" id="MathJax-Span-23" style="padding-left: 0.307em; padding-right: 0.12em;"><span style="display: inline-block; position: relative; width: 1.713em; height: 0px;"><span style="position: absolute; top: -3.08em; left: 50%; clip: rect(1.879em 1000em 2.836em -0.512em); margin-left: -0.796em;"><span class="mrow" id="MathJax-Span-24"><span class="msubsup" id="MathJax-Span-25"><span style="display: inline-block; position: relative; width: 0.814em; height: 0px;"><span style="position: absolute; top: -2.534em; left: 0em; clip: rect(2.054em 1000em 2.719em -0.512em);"><span class="mi" id="MathJax-Span-26" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">σ</span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; top: -2.428em; left: 0.407em;"><span class="texatom" id="MathJax-Span-27"><span class="mrow" id="MathJax-Span-28"><span class="mi" id="MathJax-Span-29" style="font-size: 50%; font-family: STIXGeneral; font-style: italic;">s</span><span class="mi" id="MathJax-Span-30" style="font-size: 50%; font-family: STIXGeneral; font-style: italic;">t</span></span></span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span><span class="mo" id="MathJax-Span-31" style="font-size: 70.7%; font-family: STIXGeneral;">(</span><span class="mi" id="MathJax-Span-32" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">v</span><span class="mo" id="MathJax-Span-33" style="font-size: 70.7%; font-family: STIXGeneral;">)</span></span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; top: -2.189em; left: 50%; clip: rect(2.
054em 1000em 2.823em -0.512em); margin-left: -0.413em;"><span class="msubsup" id="MathJax-Span-34"><span style="display: inline-block; position: relative; width: 0.814em; height: 0px;"><span style="position: absolute; top: -2.534em; left: 0em; clip: rect(2.054em 1000em 2.719em -0.512em);"><span class="mi" id="MathJax-Span-35" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">σ</span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; top: -2.428em; left: 0.407em;"><span class="texatom" id="MathJax-Span-36"><span class="mrow" id="MathJax-Span-37"><span class="mi" id="MathJax-Span-38" style="font-size: 50%; font-family: STIXGeneral; font-style: italic;">s</span><span class="mi" id="MathJax-Span-39" style="font-size: 50%; font-family: STIXGeneral; font-style: italic;">t</span></span></span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; top: -2.754em; left: 0em; clip: rect(2.283em 1000em 2.711em -0.531em);"><span style="border-left: 1.713em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.074em; vertical-align: 0em;"></span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 1.569em; vertical-align: -0.485em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-2">BC(v) = \sum_{s \neq v \neq t}
\frac{\sigma_{st}(v)}{\sigma_{st}}</script>, gdzie <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-40"><span style="display: inline-block; position: relative; width: 1.12em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.928em 1000em 2.868em -0.503em);"><span class="mrow" id="MathJax-Span-41"><span class="msubsup" id="MathJax-Span-42"><span style="display: inline-block; position: relative; width: 1.093em; height: 0px;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.928em 1000em 2.72em -0.503em);"><span class="mi" id="MathJax-Span-43" style="font-family: STIXGeneral; font-style: italic;">σ</span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span><span style="position: absolute; top: -2.382em; left: 0.547em;"><span class="texatom" id="MathJax-Span-44"><span class="mrow" id="MathJax-Span-45"><span class="mi" id="MathJax-Span-46" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">s</span><span class="mi" id="MathJax-Span-47" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">t</span></span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.623em; vertical-align: -0.169em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-3">\sigma_{st}</script> jest liczbą
najkrótszych ścieżek między węzłami <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-48"><span style="display: inline-block; position: relative; width: 0.413em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.914em 1000em 2.722em -0.514em);"><span class="mrow" id="MathJax-Span-49"><span class="mi" id="MathJax-Span-50" style="font-family: STIXGeneral; font-style: italic;">s</span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.482em; vertical-align: -0.014em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-4">s</script> i <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-51"><span style="display: inline-block; position: relative; width: 0.295em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.81em 1000em 2.72em -0.492em);"><span class="mrow" id="MathJax-Span-52"><span class="mi" id="MathJax-Span-53" style="font-family: STIXGeneral; font-style: italic;">t</span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.591em; vertical-align: -0.012em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-5">t</script>; a <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-54"><span style="display: inline-block; position: relative; width: 2.181em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.68em 1000em 2.886em -0.503em);"><span class="mrow" id="MathJax-Span-55"><span class="msubsup" id="MathJax-Span-56"><span style="display: inline-block; position: relative; width: 1.093em; height: 0px;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.928em 1000em 2.72em -0.503em);"><span class="mi" id="MathJax-Span-57" style="font-family: STIXGeneral; font-style: italic;">σ</span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span><span style="position: absolute; top: -2.382em; left: 0.547em;"><span class="texatom" id="MathJax-Span-58"><span class="mrow" id="MathJax-Span-59"><span class="mi" id="MathJax-Span-60" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">s</span><span class="mi" id="MathJax-Span-61" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">t</span></span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span></span><span class="mo" id="MathJax-Span-62" style="font-family: STIXGeneral;">(</span><span class="mi" id="MathJax-Span-63" style="font-family: STIXGeneral; font-style: italic;">v</span><span class="mo" id="MathJax-Span-64" style="font-family: STIXGeneral;">)</span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.905em; vertical-align: -0.188em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-6">\sigma_{st}(v)</script> jest
liczbą najkrótszych ścieżek między węzłami <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-65"><span style="display: inline-block; position: relative; width: 0.413em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.914em 1000em 2.722em -0.514em);"><span class="mrow" id="MathJax-Span-66"><span class="mi" id="MathJax-Span-67" style="font-family: STIXGeneral; font-style: italic;">s</span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.482em; vertical-align: -0.014em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-7">s</script> i <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-68"><span style="display: inline-block; position: relative; width: 0.295em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.534em; left: 0em; clip: rect(1.811em 1000em 2.722em -0.493em);"><span class="mrow" id="MathJax-Span-69"><span class="mi" id="MathJax-Span-70" style="font-family: STIXGeneral; font-style: italic;">t</span></span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.59em; vertical-align: -0.012em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-8">t</script> przechodzącymi
przez <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-71"><span style="display: inline-block; position: relative; width: 0.472em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.915em 1000em 2.727em -0.51em);"><span class="mrow" id="MathJax-Span-72"><span class="mi" id="MathJax-Span-73" style="font-family: STIXGeneral; font-style: italic;">v</span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.487em; vertical-align: -0.019em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-9">v</script>.
</p>

<p>
W sieciach opisanych grafami ważonymi obliczenie pośrednictwa wymaga
obliczenia najkrótszych ścieżek między wszystkimi parami węzłów
(np. algorytmem Floyda-Warshalla), z kosztem <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-74"><span style="display: inline-block; position: relative; width: 3.068em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -3.183em; left: 0em; clip: rect(2.125em 1000em 3.537em -0.471em);"><span class="mrow" id="MathJax-Span-75"><span class="mi" id="MathJax-Span-76" style="font-family: STIXGeneral; font-style: italic;">O</span><span class="mo" id="MathJax-Span-77" style="font-family: STIXGeneral;">(</span><span class="texatom" id="MathJax-Span-78"><span class="mrow" id="MathJax-Span-79"><span class="mo" id="MathJax-Span-80" style="vertical-align: 0.676em;"><span style="display: inline-block; position: relative; width: 0.266em; height: 0px;"><span style="position: absolute; font-family: STIXGeneral; top: -1.844em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; font-family: STIXGeneral; top: -1.871em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span></span></span><span class="mi" id="MathJax-Span-81" style="font-family: STIXGeneral; font-style: italic;">V</span><span class="msubsup" id="MathJax-Span-82"><span style="display: inline-block; position: relative; width: 0.799em; height: 0px;"><span style="position: absolute; top: -3.183em; left: 0em; clip: rect(2.33em 1000em 3.536em -0.531em);"><span class="texatom" id="MathJax-Span-83"><span class="mrow" id="MathJax-Span-84"><span class="mo" id="MathJax-Span-85" style="vertical-align: 0.676em;"><span style="display: inline-block; position: relative; width: 0.266em; height: 0px;"><span style="position: absolute; font-family: STIXGeneral; top: -1.844em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; font-family: STIXGeneral; top: -1.871em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.183em;"></span></span><span style="position: absolute; top: -2.937em; left: 0.37em;"><span class="mn" id="MathJax-Span-86" style="font-size: 70.7%; font-family: STIXGeneral;">3</span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span><span class="mo" id="MathJax-Span-87" style="font-family: STIXGeneral;">)</span></span><span style="display: inline-block; width: 0px; height: 3.183em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 1.121em; vertical-align: -0.188em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-10">O(|V|^3)</script>. Na szczęście,
wiele ciekawych sieci nie jest ważonych—do takich sieci można
stosować algorytm Brandesa (<span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-88"><span style="display: inline-block; position: relative; width: 3.658em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -3.183em; left: 0em; clip: rect(2.33em 1000em 3.537em -0.471em);"><span class="mrow" id="MathJax-Span-89"><span class="mi" id="MathJax-Span-90" style="font-family: STIXGeneral; font-style: italic;">O</span><span class="mo" id="MathJax-Span-91" style="font-family: STIXGeneral;">(</span><span class="texatom" id="MathJax-Span-92"><span class="mrow" id="MathJax-Span-93"><span class="mo" id="MathJax-Span-94" style="vertical-align: 0.676em;"><span style="display: inline-block; position: relative; width: 0.266em; height: 0px;"><span style="position: absolute; font-family: STIXGeneral; top: -1.844em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; font-family: STIXGeneral; top: -1.871em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span></span></span><span class="mi" id="MathJax-Span-95" style="font-family: STIXGeneral; font-style: italic;">V</span><span class="texatom" id="MathJax-Span-96"><span class="mrow" id="MathJax-Span-97"><span class="mo" id="MathJax-Span-98" style="vertical-align: 0.676em;"><span style="display: inline-block; position: relative; width: 0.266em; height: 0px;"><span style="position: absolute; font-family: STIXGeneral; top: -1.844em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; font-family: STIXGeneral; top: -1.871em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span></span></span><span class="texatom" id="MathJax-Span-99"><span class="mrow" id="MathJax-Span-100"><span class="mo" id="MathJax-Span-101" style="vertical-align: 0.676em;"><span style="display: inline-block; position: relative; width: 0.266em; height: 0px;"><span style="position: absolute; font-family: STIXGeneral; top: -1.844em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; font-family: STIXGeneral; top: -1.871em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span></span></span><span class="mi" id="MathJax-Span-102" style="font-family: STIXGeneral; font-style: italic;">E</span><span class="texatom" id="MathJax-Span-103"><span class="mrow" id="MathJax-Span-104"><span class="mo" id="MathJax-Span-105" style="vertical-align: 0.676em;"><span style="display: inline-block; position: relative; width: 0.266em; height: 0px;"><span style="position: absolute; font-family: STIXGeneral; top: -1.844em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; font-family: STIXGeneral; top: -1.871em; left: 0em;">∣<span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span></span></span><span class="mo" id="MathJax-Span-106" style="font-family: STIXGeneral;">)</span></span><span style="display: inline-block; width: 0px; height: 3.183em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.904em; vertical-align: -0.188em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-11">O(|V||E|)</script>).
</p>

<p>
Prosimy o zaimplementowanie algorytmu Brandesa w C++ 14.
</p>
</div>
</div>

<div id="outline-container-sec-2" class="outline-2">
<h2 id="sec-2"><span class="section-number-2">2</span> Algorytm Brandesa</h2>
<div class="outline-text-2" id="text-2">
<p>
Pseudokod współbieżnego algorytmu Brandesa przedstawiamy poniżej; dalsze wyjaśnienia w artykułach z bibliografii. Algorytm używa następujących zmiennych:
<code>BC</code>: pośrednictwo; <code>V</code>: wierzchołki; <code>d[w]</code> odległość do wierzchołka w; <code>sigma[w]</code> liczba najkrótszych ścieżek do wierzchołka w; <code>P[w]</code> poprzednicy wierzchołka <code>w</code> na wszystkich najkrótszych ścieżkach; <code>delta[v]</code> wartość pośrednictwa dla <code>v</code> w ścieżkach startujących z <code>s</code>.
</p>

<div class="org-src-container">

<pre class="src src-c"><span class="org-keyword">for</span> each v : V <br>
&nbsp;&nbsp; BC[v] = 0;<br>
<span class="org-keyword">for</span> each s : V { <span class="org-comment-delimiter">// </span><span class="org-comment">in parallel</span><br>
&nbsp;&nbsp; S = stack();<br>
&nbsp;&nbsp; <span class="org-keyword">for</span> all w in V {<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;P[w] = list();<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sigma[w] = 0;<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;d[w] = -1;<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;delta[v] = 0;<br>
&nbsp;&nbsp; }<br>
&nbsp;&nbsp; sigma[s] = 1;<br>
&nbsp;&nbsp; d[s] = 0;<br>
&nbsp;&nbsp; Q = queue(); <span class="org-comment-delimiter">// </span><span class="org-comment">FIFO</span><br>
&nbsp;&nbsp; Q.push_back(s);<br>
&nbsp;&nbsp; <span class="org-keyword">while</span> (<span class="org-negation-char">!</span>Q.empty()) {<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;v = Q.pop_front();<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;S.push(v);<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="org-keyword">for</span> each neighbor w of v {<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="org-keyword">if</span> d[w] &lt; 0 {<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Q.push_back(w);<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;d[w] = d[v] + 1;<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; }<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="org-keyword">if</span> (d[w] == d[v] + 1) {<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sigma[w] += sigma[v];<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;P[w].append(v);<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; }<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;}<br>
&nbsp;&nbsp; }<br>
&nbsp;&nbsp; <span class="org-keyword">while</span> (<span class="org-negation-char">!</span>S.empty()) {<br>
&nbsp;&nbsp;&nbsp;&nbsp; w = S.pop();<br>
&nbsp;&nbsp;&nbsp;&nbsp; <span class="org-keyword">for</span> each v in P[w]<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;delta[v] += (sigma[v] / sigma[w])(1 + delta[w]);<br>
&nbsp;&nbsp;&nbsp;&nbsp; <span class="org-keyword">if</span> (w != s)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;BC[w] += delta[w];<br>
&nbsp;&nbsp; }<br>
}<br>
<p></p>
</pre>
</div>
</div>
</div>


<div id="outline-container-sec-3" class="outline-2">
<h2 id="sec-3"><span class="section-number-2">3</span> Wejście i wyjście</h2>
<div class="outline-text-2" id="text-3">
<p>
Programy będą testowane automatycznie. Prosimy o ścisłe przestrzeganie
podanych poniżej: formatowania nazw plików oraz wejścia i wyjścia
programów.
</p>
</div>

<div id="outline-container-sec-3-1" class="outline-3">
<h3 id="sec-3-1"><span class="section-number-3">3.1</span> Sposób uruchomienia programu</h3>
<div class="outline-text-3" id="text-3-1">
<p>
<code>unzip xx123456; cd xx123456; mkdir build; cd build; cmake ..; make</code>
</p>

<p>
<code>./brandes liczba-wątków plik-wejsciowy plik-wyjsciowy</code>
</p>
</div>
</div>

<div id="outline-container-sec-3-2" class="outline-3">
<h3 id="sec-3-2"><span class="section-number-3">3.2</span> Wejście</h3>
<div class="outline-text-3" id="text-3-2">
<p>
Plik wejściowy opisuje krawędzie w grafie w formacie <code>&lt;wierzchołek1&gt; &lt;wierzchołek2&gt;</code>. Graf jest skierowany. Krawędzie są posortowane po
pierwszym wierzchołku i po drugim wierzchołku. Wierzchołki nie muszą
być numerowane sekwencyjnie (tzn. poprawny jest graf bez
np. wierzchołka o numerze 0). Graf nie musi być spójny. 
</p>

<p>
Program powinien przeprowadzać obliczenia używając <code>liczba-wątków</code>
wątków. W testach wydajnościowych będzie zdecydowanie więcej
wierzchołków niż wątków, ale musisz zapewnić poprawne działanie
programu dla <code>liczba-wątków &lt; 100</code> niezależnie od liczby wierzchołków.
</p>


<p>
Przykładowe wejście:
</p>
<pre class="example">0 2
2 0
2 3
2 4
3 2
3 5
3 6
</pre>
</div>
</div>

<div id="outline-container-sec-3-3" class="outline-3">
<h3 id="sec-3-3"><span class="section-number-3">3.3</span> Wyjście</h3>
<div class="outline-text-3" id="text-3-3">
<p>
Plik wyjściowy składa się z dokładnie tylu linii, ile jest
wierzchołków z co najmniej jedną krawędzią wychodzącą. <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-107"><span style="display: inline-block; position: relative; width: 0.295em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.702em 1000em 2.72em -0.481em);"><span class="mrow" id="MathJax-Span-108"><span class="mi" id="MathJax-Span-109" style="font-family: STIXGeneral; font-style: italic;">i</span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.705em; vertical-align: -0.012em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-12">i</script>-ta linia
zawiera dwie liczby: liczbę całkowitą — identyfikator
wierzchołka, oraz liczbę rzeczywistą — pośrednictwo tego
wierzchołka. Wyjście musi być posortowane po pierwszej kolumnie.
</p>

<p>
Standardowe wyjście programu będzie ignorowane.
</p>

<p>
Plik wyjściowy odpowiadający przykładowemu plikowi wejściowemu:
</p>

<pre class="example">0 0
2 6
3 4
</pre>
</div>
</div>
</div>


<div id="outline-container-sec-4" class="outline-2">
<h2 id="sec-4"><span class="section-number-2">4</span> Forma oddania zadania</h2>
<div class="outline-text-2" id="text-4">
<p>
Prosimy o oddanie pojedynczego pliku <code>.zip</code> zawierającego pojedynczy katalog odpowiadający loginowi (<code>ab123456</code>), a w nim następujące pliki:
</p>
<ul class="org-ul">
<li>kod źródłowy;
</li>
<li>CMakeLists.txt;
</li>
<li>raport w pliku raport.pdf
</li>
</ul>

<p>
W raporcie opisz twoją implementację i zastosowane optymalizacje. Dla
grafu <code>wiki-vote-sort.txt</code> zmierz speed-up, czyli zmierz czas
działania twojego programu <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-110"><span style="display: inline-block; position: relative; width: 1.239em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.534em; left: 0em; clip: rect(1.916em 1000em 2.958em -0.606em);"><span class="mrow" id="MathJax-Span-111"><span class="msubsup" id="MathJax-Span-112"><span style="display: inline-block; position: relative; width: 1.212em; height: 0px;"><span style="position: absolute; top: -2.534em; left: 0em; clip: rect(1.916em 1000em 2.916em -0.606em);"><span class="mi" id="MathJax-Span-113" style="font-family: STIXGeneral; font-style: italic;">p</span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span><span style="position: absolute; top: -2.294em; left: 0.606em;"><span class="mi" id="MathJax-Span-114" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">m</span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 2.534em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.729em; vertical-align: -0.261em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-13">p_m</script> dla <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-115"><span style="display: inline-block; position: relative; width: 6.071em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.68em 1000em 2.85em -0.518em);"><span class="mrow" id="MathJax-Span-116"><span class="mi" id="MathJax-Span-117" style="font-family: STIXGeneral; font-style: italic;">m</span><span class="mo" id="MathJax-Span-118" style="font-family: STIXGeneral; padding-left: 0.313em;">=</span><span class="mn" id="MathJax-Span-119" style="font-family: STIXGeneral; padding-left: 0.313em;">1</span><span class="mo" id="MathJax-Span-120" style="font-family: STIXGeneral;">,</span><span class="mn" id="MathJax-Span-121" style="font-family: STIXGeneral; padding-left: 0.188em;">2</span><span class="mo" id="MathJax-Span-122" style="font-family: STIXGeneral;">,</span><span class="mo" id="MathJax-Span-123" style="font-family: STIXGeneral; padding-left: 0.188em;">.</span><span class="mo" id="MathJax-Span-124" style="font-family: STIXGeneral; padding-left: 0.188em;">.</span><span class="mo" id="MathJax-Span-125" style="font-family: STIXGeneral; padding-left: 0.188em;">.</span><span class="mo" id="MathJax-Span-126" style="font-family: STIXGeneral; padding-left: 0.188em;">,</span><span class="mn" id="MathJax-Span-127" style="font-family: STIXGeneral; padding-left: 0.188em;">8</span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.866em; vertical-align: -0.15em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-14">m=1, 2, ..., 8</script> wątków, a
następnie oblicz przyspieszenie: <span class="MathJax_Preview" style="display: none;"></span><span class="MathJax" role="textbox" aria-readonly="true"><nobr><span class="math" id="MathJax-Span-128"><span style="display: inline-block; position: relative; width: 4.951em; height: 0px; font-size: 106%;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.68em 1000em 2.956em -0.513em);"><span class="mrow" id="MathJax-Span-129"><span class="msubsup" id="MathJax-Span-130"><span style="display: inline-block; position: relative; width: 1.152em; height: 0px;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.689em 1000em 2.727em -0.513em);"><span class="mi" id="MathJax-Span-131" style="font-family: STIXGeneral; font-style: italic;">S</span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span><span style="position: absolute; top: -2.382em; left: 0.547em;"><span class="mi" id="MathJax-Span-132" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">m</span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span></span><span class="mo" id="MathJax-Span-133" style="font-family: STIXGeneral; padding-left: 0.313em;">=</span><span class="msubsup" id="MathJax-Span-134" style="padding-left: 0.313em;"><span style="display: inline-block; position: relative; width: 1.034em; height: 0px;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.915em 1000em 2.914em -0.605em);"><span class="mi" id="MathJax-Span-135" style="font-family: STIXGeneral; font-style: italic;">p</span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span><span style="position: absolute; top: -2.292em; left: 0.605em;"><span class="mn" id="MathJax-Span-136" style="font-size: 70.7%; font-family: STIXGeneral;">1</span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span></span><span class="texatom" id="MathJax-Span-137"><span class="mrow" id="MathJax-Span-138"><span class="mo" id="MathJax-Span-139" style="font-family: STIXGeneral;">/</span></span></span><span class="msubsup" id="MathJax-Span-140"><span style="display: inline-block; position: relative; width: 1.211em; height: 0px;"><span style="position: absolute; top: -2.532em; left: 0em; clip: rect(1.915em 1000em 2.914em -0.605em);"><span class="mi" id="MathJax-Span-141" style="font-family: STIXGeneral; font-style: italic;">p</span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span><span style="position: absolute; top: -2.292em; left: 0.605em;"><span class="mi" id="MathJax-Span-142" style="font-size: 70.7%; font-family: STIXGeneral; font-style: italic;">m</span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 2.532em;"></span></span></span><span style="border-left: 0em solid; display: inline-block; overflow: hidden; width: 0px; height: 0.978em; vertical-align: -0.262em;"></span></span></nobr></span><script type="math/tex" id="MathJax-Element-15">S_m=p_1 / p_m</script>.
</p>

<p>
CMakeLists powinien kompilować program do jak najszybszej wersji:
pamiętaj o włączeniu optymalizacji kompilatora i wyłączeniu debugowania.
</p>
</div>
</div>

<div id="outline-container-sec-5" class="outline-2">
<h2 id="sec-5"><span class="section-number-2">5</span> Kryteria oceny</h2>
<div class="outline-text-2" id="text-5">
<p>
Będziemy oceniać poprawność, wydajność oraz raport. 
</p>

<p>
Poprawność będziemy oceniać porównując wyniki implementacji z naszym
rozwiązaniem wzorcowym; nasze testy biorą pod uwagę błędy
zaokrągleń. Udostępniamy dwa przykładowe testy, ale nie wyczerpują one
wszystkich przypadków brzegowych — pamiętaj o dodaniu własnych
testów. Nie będziemy oceniać wydajności programów niepoprawnych.
</p>

<p>
Szybkość działania będziemy oceniać porównując czas działania waszych
implementacji na sieciach "społecznych", tzn. nie
Erdős'a–Rényi'ego. Pojedyncze wykonanie programu będziemy ograniczać
do kilku (max. 5) minut.
</p>
</div>
</div>


<div id="outline-container-sec-6" class="outline-2">
<h2 id="sec-6"><span class="section-number-2">6</span> Dodatkowe materiały</h2>
<div class="outline-text-2" id="text-6">
<p>
Prosimy o niekorzystanie z kodów źródłowych - gotowych implementacji
algorytmu Brandesa. Zalecamy natomiast przeczytanie następujących prac:
</p>

<ul class="org-ul">
<li>A faster algorithm for betweenness centrality, U Brandes, Journal of
Mathematical Sociology, 2001 Taylor &amp; Francis, <a href="https://kops.uni-konstanz.de/bitstream/handle/123456789/5739/algorithm.pdf?sequence=1">https://kops.uni-konstanz.de/bitstream/handle/123456789/5739/algorithm.pdf?sequence=1</a>
</li>
<li>Parallel algorithms for evaluating centrality indices in real-world
networks DA Bader, K Madduri, <a href="https://smartech.gatech.edu/bitstream/handle/1853/14428/GT-CSE-06-13.pdf">https://smartech.gatech.edu/bitstream/handle/1853/14428/GT-CSE-06-13.pdf</a>
</li>
<li>Dane do testów wydajnościowych można wziąć z <a href="http://snap.stanford.edu/data/">http://snap.stanford.edu/data/</a>.
</li>
</ul>
</div>
</div>


<div id="outline-container-sec-7" class="outline-2">
<h2 id="sec-7"><span class="section-number-2">7</span> FAQ</h2>
<div class="outline-text-2" id="text-7">
<ul class="org-ul">
<li>Dokładność obliczeń
</li>
</ul>

<p>
Proszę używać liczb zmiennoprzecinkowych o podwójnej precyzji.
</p>

<ul class="org-ul">
<li>Czy możemy założyc, że wejście jest poprawne?
</li>
</ul>

<p>
Tak. Czyli program może zachowywać się nieprzewidywalnie dla wejść
niespełniających specyfikacji.
</p>

<ul class="org-ul">
<li>Czy możemy stworzyć dodatkowe, pomocnicze wątki?
</li>
</ul>

<p>
Liczba wątków "liczących", tzn. wykonujących algorytm Bradesa, musi
być równa <code>liczba-watkow</code>, pierwszemu parametrowi wywołania programu. Można stworzyć
dodatkowe, pomocnicze wątki, ale co najwyżej O(<code>liczba-watkow</code>). Dodatkowe wątki
mogą wykonać co najwyżej O(1) operacji.
</p>

<ul class="org-ul">
<li>Czy można korzystać z dodatkowych bibliotek?
</li>
</ul>

<p>
Nie. W rozwiązaniu można korzystać tylko ze standardowej biblioteki
CPP14. W szczególności nie można korzystać z Boosta.
</p>

<ul class="org-ul">
<li>Co można założyć o identyfikatorach węzłów grafu?
</li>
</ul>

<p>
Zmieszczą się do zwykłego <code>int</code>.
</p>
</div>
</div>
</div>
<div id="postamble" class="status">
<p class="date">Data: 2016/12/26</p>
<p class="author">Autor: Krzysztof Rządca</p>
<p class="date">Created: 2016-12-26 Mon 21:06</p>
<p class="creator"><a href="http://www.gnu.org/software/emacs/">Emacs</a> 25.1.1 (<a href="http://orgmode.org/">Org</a> mode 8.2.10)</p>
<p class="validation"><a href="http://validator.w3.org/check?uri=referer">Validate</a></p>
</div>


<div style="position: absolute; visibility: hidden; top: 0px; left: 0px; font-size: 40px; font-weight: normal; font-style: normal; font-family: STIXSizeOneSym, sans-serif;">ˆ</div></body></html>

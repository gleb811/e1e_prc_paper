#!/bin/tcsh -f
rm ad_hoc_first.pdf
rm ad_hoc_first.aux
#rm text/aliases/*.aux text/conclusions/*.aux text/event_selection/*.aux text/motivation/*.aux text/setups/*.aux text/title_page/*.aux text/toc/*.aux text/appendix/*.aux text/rad_corr/rad_corr.aux text/sys_err/*.aux text/bin_corr/*.aux text/cross_sction/*.aux
pdflatex ad_hoc_first
pdflatex ad_hoc_first
acroread ad_hoc_first.pdf >& /dev/null
rm ad_hoc_first.log ad_hoc_first.aux ad_hoc_first.toc ad_hoc_first.out 

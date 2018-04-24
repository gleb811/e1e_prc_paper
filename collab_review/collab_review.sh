#!/bin/tcsh -f
rm collab_review.pdf
rm collab_review.aux
#rm text/aliases/*.aux text/conclusions/*.aux text/event_selection/*.aux text/motivation/*.aux text/setups/*.aux text/title_page/*.aux text/toc/*.aux text/appendix/*.aux text/rad_corr/rad_corr.aux text/sys_err/*.aux text/bin_corr/*.aux text/cross_sction/*.aux
pdflatex collab_review
pdflatex collab_review
acroread collab_review.pdf >& /dev/null
rm collab_review.log collab_review.aux collab_review.toc collab_review.out 

#!/bin/bash
cd /work/users/s/e/seyoun/chon_page
Rscript -e "shiny::runApp(host='0.0.0.0', port=8888)"

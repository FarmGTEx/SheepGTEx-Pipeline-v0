jsub -q normal -n 1 -R "span[hosts=1]" -J sheep_combine -e log/sheep_combine.%J.log -o log/sheep_combine.%J.log "python combine.py"

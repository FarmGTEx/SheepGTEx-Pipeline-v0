#!/bin/bash
awk '{print $1}' chrAuto.filtered.keep.fam | csvtk join -t -H -f '1;1' - ../group.list > index.list
paste index.list chrAuto.filtered.keep.pruned.2.Q > chrAuto.filtered.keep.pruned.merged.Q
for k in {3..10} ; do paste chrAuto.filtered.keep.pruned.merged.Q chrAuto.filtered.keep.pruned.${k}.Q > chrAuto.filtered.keep.pruned.merged.Q.tmp ; mv chrAuto.filtered.keep.pruned.merged.Q.tmp chrAuto.filtered.keep.pruned.merged.Q ; done
sed -i 's/ /\t/g' chrAuto.filtered.keep.pruned.merged.Q

PIPELINE_TOOLS=$(echo tools/* | tr ' ' '\n' | xargs -I'{}' bash -c 'td=$(readlink -f $1)/current; bd=$1/.bindir; [[ -f $bd ]] && echo $td/$(cat $bd) || echo $td' -- {} | tr '\n' ':' | sed -re 's/.$//')
export PATH=$PIPELINE_TOOLS:$PATH

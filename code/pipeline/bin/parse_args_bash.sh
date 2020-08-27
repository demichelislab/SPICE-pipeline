i=$(($# + 1))
[[ -v longoptspec[@] ]] || declare -A longoptspec
[[ -v optspec ]] || optspec=""
if (( ${#optspec} == 0 && ${#longoptspec[@]} == 0 ))
then
  set_args () {
    false
  }
fi
longoptspec[help]=0
optspec=":${optspec}h-:"
while getopts "${optspec}" opt; do
while true; do
    if [[ "${opt}" != "-" ]]
    then
      if set_args "${opt}" "${OPTARG}"
      then
        break;
      fi
    fi
    case "${opt}" in
        -) #OPTARG is name-of-long-option or name-of-long-option=value
            if [[ ${OPTARG} =~ .*=.* ]] # with this --key=value format only one argument is possible
            then
                opt=${OPTARG/=*/}
                ((${#opt} <= 1)) && {
                    echo "Syntax error: Invalid long option '$opt'" >&2
                    exit 2
                }
                if (($((longoptspec[$opt])) != 1))
                then
                    echo "Syntax error: Option '$opt' does not support this syntax." >&2
                    exit 2
                fi
                OPTARG=${OPTARG#*=}
            else #with this --key value1 value2 format multiple arguments are possible
                opt="$OPTARG"
                ((${#opt} <= 1)) && {
                    echo "Syntax error: Invalid long option '$opt'" >&2
                    exit 2
                }
                OPTARG=(${@:OPTIND:$((longoptspec[$opt]))})
                ((OPTIND+=longoptspec[$opt]))
                ((OPTIND > i)) && {
                    echo "Syntax error: Not all required arguments for option '$opt' are given." >&2
                    exit 3
                }
            fi
            continue #now that opt/OPTARG are set we can process them as
            # if getopts would've given us long options
            ;;
        h|help)
            usage
            exit 1
            ;;
        ?)
            echo "Syntax error: Unknown short option '$OPTARG'" >&2
            exit 2
            ;;
        *)
            echo "Syntax error: Unknown long option '$opt'" >&2
            exit 2
            ;;
    esac
break; done
done
shift $((OPTIND -1))
unset i optspec longoptspec
if [[ $# -lt ${required_args} ]]
then
  usage
  exit 1
fi


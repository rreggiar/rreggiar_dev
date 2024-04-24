#!/usr/bin/env bash
#   Installs IgBlast for use with TCR package for split-pipe pipeline
#
#   Copyright (c) 2023 Parse Biosciences
#
#   See company page for background information, usage instructions and license.
#   https://www.parsebiosciences.com/
#

set -o errexit
set -o nounset


VERSION="Version 0.4; RTK,AK 2023-11-11"

# For version, save program name from first var, first pass
PROGNAME=$0
version() { echo "$PROGNAME $VERSION"; }

# Output top level
igb_top="IgBlast"

# Log file name; Full path below
PROC_LOG_NAME="install_IgBlast.log"


# BLAST Target url files
FTP_URL_LATEST="https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/"
FTP_URL_STABLE="https://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.20.0/"
IDX_FILE="index.html"

# IMGT database collection and download command (wget)
imgt_targs="\
https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRAV.fasta \
https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRAJ.fasta \
https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBV.fasta \
https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBD.fasta \
https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBJ.fasta \
https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP \
"
# May not work with simple wget; Need extra arg
imgt_wget="wget"
#imgt_wget="wget --no-check-certificate"


main ()
{   
    if [[ $# -lt 4 ]]; then
        if [[ $# -lt 1 ]]; then
            echo "Usage: <path_above_IgBlast> [<latest>]"
            echo " "
            echo "If second arg, install from 'LATEST' ftp path"
            exit 1
        fi

        # Target source URL and version indication
        if [[ $2 == 'latest' ]]; then
            FTP_URL="$FTP_URL_LATEST"
        else
            FTP_URL="$FTP_URL_STABLE"
        fi
        TARG_VERSION=$(echo $FTP_URL | tr "/" " " | awk '{print $NF}')
        
        top_dir=$1
        if [[ ! -d $top_dir ]]; then 
            echo "Bad top level dir given: $top_dir"
            exit 1
        fi
        # Get full path for top_dir (in case given relative)
        top_dir=$(cd "$top_dir" || exit; pwd)

        # Save processing outputs to log file
        PROC_LOG="${top_dir}/${PROC_LOG_NAME}"
        [[ -e $PROC_LOG ]] && rm "$PROC_LOG"
        {
            echo "Install for IgBlast"
            date
            echo "-----------------"
            echo " "                    
        } >> "$PROC_LOG"

        # Save starting location and jump to top dir
        start_top=$(pwd)
        cd "$top_dir" || exit

        echo "# Saving details to: $PROC_LOG"
        echo "# Will install $igb_top under $top_dir" | tee -a "$PROC_LOG"
        echo "# Target version $TARG_VERSION" | tee -a "$PROC_LOG"

        n_start_index=$(find . -name "index.html" | wc | awk '{print $1}')
        echo "# N start index $n_start_index" | tee -a "$PROC_LOG"

        ####  Install igBLAST package
        # Figure out lastest version download command
        # First get index to find latest tar version name from that
        echo "# Getting index for $FTP_URL" | tee -a "$PROC_LOG"
        call="wget $FTP_URL"
        {
            echo "$call"
            eval "$call"
            echo "-----------------"
            echo " "
        } &>> "$PROC_LOG"

        if [[ ! -f $IDX_FILE ]]; then
            echo "Failed to get $IDX_FILE from $FTP_URL" | tee -a "$PROC_LOG"
            exit 1
        fi

        # Get target ftp and filename from html index
        #ftp_tar_file=$(grep x64-linux.tar $IDX_FILE | grep -v md5 | sed 's/[<>=]/ /g' | awk '{print $(NF-4)}' | sed 's/"//g')
        #tar_file=$(grep x64-linux.tar $IDX_FILE | grep -v md5 | sed 's/[<>=]/ /g' | awk '{print $(NF-3)}')
        tar_file=$(grep x64-linux.tar $IDX_FILE | grep -v md5 | sed 's/[<>=]/ /g' | awk '{print $(NF-4)}')
        ftp_tar_file=$FTP_URL$tar_file
        echo "# Getting tarball $ftp_tar_file" | tee -a "$PROC_LOG"
        call="wget $ftp_tar_file"
        {
            echo "$call"
            eval "$call"
            echo "-----------------"
            echo " "
        } &>> "$PROC_LOG"

        if [[ -f $tar_file ]]; then
            echo "# Got $tar_file from $FTP_URL" | tee -a "$PROC_LOG"
        else
            echo "Failed to get $tar_file from $FTP_URL" | tee -a "$PROC_LOG"
            exit 1
        fi

        # Remove 'index.html' if this is new
        if [[ -f $IDX_FILE  ]] && [[ $n_start_index -eq 0 ]]; then
            rm $IDX_FILE
        fi

        # Untar then get dir name
        echo "# Extracting $tar_file" | tee -a "$PROC_LOG"
        call="tar -xf $tar_file"
        {
            echo "$call"
            eval "$call"
            echo "-----------------"
            echo " "
        } &>> "$PROC_LOG"

        raw_igb_top=$(find . -maxdepth 1 -type d -name "*igblast*")
        if [[ -d $raw_igb_top ]]; then
            echo "# Extracted $raw_igb_top" | tee -a "$PROC_LOG"
        else
            echo "Problem extracting from $tar_file" | tee -a "$PROC_LOG"
            exit 1
        fi

        # Done with tar
        rm "$tar_file"

        echo "# Moving raw IgBlast $raw_igb_top >--> $igb_top" | tee -a "$PROC_LOG"
        mv "$raw_igb_top" "$igb_top"

        # Jump into top dir, make dir for fasta files and retrieve
        cd $igb_top || exit
        new_igb_path=$(pwd)

        # Make sure blast works
        call=" bin/makeblastdb -version"
        # run_blast=$( $call )
        exit_code=$?
        if [[ $exit_code -gt 0 ]]; then
            echo "Problem with call: '${call}'" | tee -a "$PROC_LOG"
            echo "Exit code $exit_code" | tee -a "$PROC_LOG"
            exit $exit_code
        fi
    else
        start_top=$(pwd)
        new_igb_path=${igb_top}
        top_dir=$1
        if [[ ! -d $top_dir ]]; then 
            echo "Bad top level dir given: $top_dir"
            exit 1
        fi
        top_dir=$(cd "$top_dir" || exit; pwd)
        # Save processing outputs to log file
        PROC_LOG="${top_dir}/${PROC_LOG_NAME}"
    fi

    ####  Install data files
    db_type=$3
    genome_name="human"
    seq_type="tcr"
    # starting point
    cd "${top_dir}/${igb_top}"

    #Check if default db exists and not empty
    raw_check=0
    if [[ -d tr_database ]]; then
        count=$(find tr_database -type f -name "${genome_name}_${seq_type}*.nsq" | wc -l)
        if [[ "$count" -eq 4 ]];then
            raw_check=1
        fi
    fi

    # Unpack and copy our files
    if [[ raw_check -eq 0 ]]; then
        echo "# Building default database"
        mkdir RAW_TR
        cd RAW_TR || exit
        if [[ "$genome_name" == "human" ]]; then
            echo "# Unzipping db to igblast directory" | tee -a "$PROC_LOG"
            cp "../../../scripts/db/${genome_name}_${seq_type}.zip" ./
            unzip "${genome_name}_${seq_type}.zip"
            mv "${genome_name}_${seq_type}/${genome_name}_${seq_type}_gl.aux" ../optional_file/
            mv ../optional_file/human_gl.aux ../optional_file/imgt_human_tcr_gl.aux
            echo "-----------------"       
        fi

        # Convert IMGT header format 
        cd ../
        mkdir tr_fasta
        bin/edit_imgt_file.pl "RAW_TR/${genome_name}_${seq_type}/${genome_name}_${seq_type}_V" > "tr_fasta/${genome_name}_${seq_type}_V"
        bin/edit_imgt_file.pl "RAW_TR/${genome_name}_${seq_type}/${genome_name}_${seq_type}_D" > "tr_fasta/${genome_name}_${seq_type}_D"
        bin/edit_imgt_file.pl "RAW_TR/${genome_name}_${seq_type}/${genome_name}_${seq_type}_J" > "tr_fasta/${genome_name}_${seq_type}_J"
        bin/edit_imgt_file.pl "RAW_TR/${genome_name}_${seq_type}/${genome_name}_${seq_type}_C" > "tr_fasta/${genome_name}_${seq_type}_C"

        echo "# Building BLAST database" | tee -a "$PROC_LOG"
        # Make the database
        mkdir tr_database
        {
            bin/makeblastdb -parse_seqids -dbtype nucl -in "tr_fasta/${genome_name}_${seq_type}_V" -out "tr_database/${genome_name}_${seq_type}_V"
            bin/makeblastdb -parse_seqids -dbtype nucl -in "tr_fasta/${genome_name}_${seq_type}_D" -out "tr_database/${genome_name}_${seq_type}_D"
            bin/makeblastdb -parse_seqids -dbtype nucl -in "tr_fasta/${genome_name}_${seq_type}_J" -out "tr_database/${genome_name}_${seq_type}_J"
            bin/makeblastdb -parse_seqids -dbtype nucl -in "tr_fasta/${genome_name}_${seq_type}_C" -out "tr_database/${genome_name}_${seq_type}_C"    
        } &>> "$PROC_LOG"
    else
        echo "Default database already exists for ${genome_name} ${seq_type}"
    fi 


    # If customer wants IMGT. Fetch and index IMGT
    if [[ "$db_type" == "imgt" ]]; then
        # starting dir
        cd "${top_dir}/${igb_top}"

        # Check if imgt db exists and not empty
        im_check=0
        if [[ -d tr_database ]]; then
            count=$(find tr_database -type f -name "imgt_${genome_name}_${seq_type}*.nsq" | wc -l)
            if [[ "$count" -eq 4 ]]; then
                im_check=1
            fi
        fi

        if [[ im_check -eq 0 ]]; then
            echo "-----------------"       
            echo "# Will add IMGT database"

            imgt_dir="IMGT_TR"
            if [[ ! -d $imgt_dir ]]; then 
                mkdir $imgt_dir
            fi
            cd $imgt_dir || exit

            # Only human for now
            if [[ "$genome_name" == "human" ]]; then

                imgt_g_dir="${genome_name}_${seq_type}"
                if [[ ! -d $imgt_g_dir ]]; then 
                    mkdir $imgt_g_dir
                fi
                cd $imgt_g_dir || exit

                echo "# Getting IMGT fasta data files" | tee -a "$PROC_LOG"
                echo "# Fetch command: $imgt_wget"
                for imgt_targ_url in $imgt_targs; do
                    echo "# Getting target $imgt_targ_url" | tee -a "$PROC_LOG"
                    # Fetch call (i.e. wget)
                    if $imgt_wget "$imgt_targ_url" &>> "$PROC_LOG"; then
                        echo "#            Got $imgt_targ_url" | tee -a "$PROC_LOG"
                    else
                        echo "$imgt_wget failed with exit code $?" | tee -a "$PROC_LOG"
                        exit 1
                    fi
                done

                echo "-----------------"
                echo "# Processing fasta files" | tee -a "$PROC_LOG"

                #Process the fasta file to get C genes
                awk -F'>' 'NF>1{f=($2 ~ /.*TR[AB]C.*Homo sapien.*/)} f' IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP > ${genome_name}_${seq_type}_C_sequences.fasta
                awk -F"|" '/TR/ {id = $2;if (pid == id) {} else if (!pid){print ">"id} else {print ">"id};pid=id} /^[acgtn.]/ {print}' ${genome_name}_${seq_type}_C_sequences.fasta> ${genome_name}_${seq_type}_C.fasta

                # Combine other segments
                {
                    cat ./*V.fasta > "${genome_name}_${seq_type}_V.fasta"
                    cat ./*D.fasta > "${genome_name}_${seq_type}_D.fasta"
                    cat ./*J.fasta > "${genome_name}_${seq_type}_J.fasta"
                } &>> "$PROC_LOG"
            else
                echo "Unknown genome $genome_name"
                exit 1
            fi

            # For now only get alpha and beta chains.
            # Note that using alpha and delta chains together will cause an error because the IMGT TRDV.fasta file
            # has TRAV sequences.
            # wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRDV.fasta
            # wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRDD.fasta
            # wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRDJ.fasta
            # wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRGV.fasta
            # wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRGJ.fasta

            # Convert IMGT files
            cd ../../
            if [[ ! -d tr_fasta ]]; then
                mkdir tr_fasta
            fi
            
            bin/edit_imgt_file.pl "IMGT_TR/${genome_name}_${seq_type}/${genome_name}_${seq_type}_V.fasta" > "tr_fasta/imgt_${genome_name}_${seq_type}_V"
            bin/edit_imgt_file.pl "IMGT_TR/${genome_name}_${seq_type}/${genome_name}_${seq_type}_D.fasta" > "tr_fasta/imgt_${genome_name}_${seq_type}_D"
            bin/edit_imgt_file.pl "IMGT_TR/${genome_name}_${seq_type}/${genome_name}_${seq_type}_J.fasta" > "tr_fasta/imgt_${genome_name}_${seq_type}_J"
            bin/edit_imgt_file.pl "IMGT_TR/${genome_name}_${seq_type}/${genome_name}_${seq_type}_C.fasta" > "tr_fasta/imgt_${genome_name}_${seq_type}_C"

            echo "# Building IMGT BLAST database" | tee -a "$PROC_LOG"
            # Make the database
            if [[ ! -d tr_database ]]; then
                mkdir tr_database
            fi
            
            {
                bin/makeblastdb -parse_seqids -dbtype nucl -in "tr_fasta/imgt_${genome_name}_${seq_type}_V" -out "tr_database/imgt_${genome_name}_${seq_type}_V"
                bin/makeblastdb -parse_seqids -dbtype nucl -in "tr_fasta/imgt_${genome_name}_${seq_type}_D" -out "tr_database/imgt_${genome_name}_${seq_type}_D"
                bin/makeblastdb -parse_seqids -dbtype nucl -in "tr_fasta/imgt_${genome_name}_${seq_type}_J" -out "tr_database/imgt_${genome_name}_${seq_type}_J"
                bin/makeblastdb -parse_seqids -dbtype nucl -in "tr_fasta/imgt_${genome_name}_${seq_type}_C" -out "tr_database/imgt_${genome_name}_${seq_type}_C"    
            } &>> "$PROC_LOG"

        else
            echo "IMGT database already exists for ${genome_name} ${seq_type}"
        fi 
    fi
    # Jump back to starting dir; 
    cd "$start_top" || exit

    echo "# IgBlast dir: $new_igb_path" | tee -a "$PROC_LOG"
    echo "# IgBlast log: $PROC_LOG"
    echo "# All done. Finished successfully" | tee -a "$PROC_LOG"

}


# Call main with all command line args
main "$@"


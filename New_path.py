
def get_transcript_scores(  # pylint: disable=too-many-locals
        s_exon_table,
        transcript_table,
        graph,
        column_order=None,
        delim='/'):
    """
    Return a DataFrame with the needed data to choose the canonical path.
    """
    data = {
        'GeneID': [],
        'TranscriptIDCluster': [],
        'TranscriptLength': [],
        'MinimumTranscriptWeightedConservation': [],
        'MinimumTranscriptFraction': [],
        'MinimumConservation': [],
        'TranscriptWeightedConservationSum': [],
        'TranscriptFractionSum': [],
        'ConservationSum': [],
        'MeanTranscriptWeightedConservation': [],
        'MeanTranscriptFraction': [],
        'MeanConservation': [],
        'TSL': [],
        'Path': []
    }

    transcript2tsl = transcript_table.loc[:, [
        "TranscriptIDCluster", "TSL"
    ]].set_index("TranscriptIDCluster").to_dict()['TSL']

    for (trx, subdf) in s_exon_table.groupby('TranscriptIDCluster'):
        n_rows = len(subdf)
        s_exon_len = [len(_get_seq(subdf.S_exon_Sequence.iloc[0]))]
        path = ["start", subdf.S_exonID.iloc[0]]
        score = [
            graph.get_edge_data(
                'start',
                subdf.S_exonID.iloc[0])['transcript_weighted_conservation']
        ]
        conservation = [
            graph.get_edge_data('start',
                                subdf.S_exonID.iloc[0])['conservation']
        ]
        transcript_fraction = [
            graph.get_edge_data('start',
                                subdf.S_exonID.iloc[0])['transcript_fraction']
        ]
        if n_rows >= 2:
            for i in range(1, n_rows):
                s_exon_1 = subdf.S_exonID.iloc[i - 1]
                s_exon_2 = subdf.S_exonID.iloc[i]
                s_exon_len.append(len(_get_seq(subdf.S_exon_Sequence.iloc[i])))
                path.append(subdf.S_exonID.iloc[i])
                score.append(
                    graph.get_edge_data(
                        s_exon_1,
                        s_exon_2)['transcript_weighted_conservation'])
                conservation.append(
                    graph.get_edge_data(s_exon_1, s_exon_2)['conservation'])
                transcript_fraction.append(
                    graph.get_edge_data(s_exon_1,
                                        s_exon_2)['transcript_fraction'])

        score.append(
            graph.get_edge_data(subdf.S_exonID.iloc[n_rows - 1],
                                'stop')['transcript_weighted_conservation'])
        conservation.append(
            graph.get_edge_data(subdf.S_exonID.iloc[n_rows - 1],
                                'stop')['conservation'])
        transcript_fraction.append(
            graph.get_edge_data(subdf.S_exonID.iloc[n_rows - 1],
                                'stop')['transcript_fraction'])
        path.append('stop')

        data['GeneID'].append(subdf.GeneID.iloc[0])
        data['TranscriptIDCluster'].append(trx)
        data['TranscriptLength'].append(sum(s_exon_len))
        data['MinimumTranscriptWeightedConservation'].append(min(score))
        data['MinimumConservation'].append(min(conservation))
        data['MinimumTranscriptFraction'].append(min(transcript_fraction))
        data['TranscriptWeightedConservationSum'].append(sum(score))
        data['ConservationSum'].append(sum(conservation))
        data['TranscriptFractionSum'].append(sum(transcript_fraction))
        data['MeanTranscriptWeightedConservation'].append(np.mean(score))
        data['MeanConservation'].append(np.mean(conservation))
        data['MeanTranscriptFraction'].append(np.mean(transcript_fraction))
        data['TSL'].append(transcript2tsl[trx])
        data['Path'].append(delim.join(path))

    data_frame = pd.DataFrame(data)
    path2ngenes = data_frame.groupby('Path').apply(
        lambda df: len(set(df.GeneID))).to_dict()
    data_frame = data_frame.assign(
        PathGeneNumber=[path2ngenes[path] for path in data_frame.Path])

    if column_order is None:
        column_order = [
            'MinimumConservation', 'MinimumTranscriptWeightedConservation',
            'MeanTranscriptWeightedConservation', 'TranscriptLength', 'TSL'
        ]

    data_frame.drop_duplicates(inplace=True)
    data_frame.sort_values(column_order, ascending=False, inplace=True)

    return data_frame.reindex(columns=[
        'GeneID', 'TranscriptIDCluster', 'TranscriptLength', 'TSL', 'Path',
        'MinimumConservation', 'MinimumTranscriptFraction',
        'MinimumTranscriptWeightedConservation', 'ConservationSum',
        'TranscriptFractionSum', 'TranscriptWeightedConservationSum',
        'MeanConservation', 'MeanTranscriptFraction',
        'MeanTranscriptWeightedConservation', 'PathGeneNumber'
    ])

def transcript_path_modification() :
    
    
    return 
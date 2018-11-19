'''
Created on Feb 21, 2018

@author: pwalters
'''
import os
import tasks
import pypeliner
import pypeliner.managed as mgd


def create_destruct_workflow(
        tumour_bam, normal_bam, raw_breakpoints, raw_library,
        breakpoints, library, reads,
        sample_id, global_config, sv_config):

    destruct_config = {}

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='destruct',
        func='destruct.workflow.create_destruct_workflow',
        args=(
            {sample_id: mgd.InputFile(tumour_bam),
             sample_id + 'N': mgd.InputFile(normal_bam)},
            mgd.OutputFile(raw_breakpoints),
            mgd.OutputFile(raw_library),
            mgd.OutputFile(reads),
            destruct_config,
            sv_config['refdata_destruct']
        )
    )

    workflow.transform(
        name='filter_annotate_breakpoints',
        ctx={
            'mem': 4,
            'pool_id': 'standard',
            'ncpus': 1,
        },
        func=tasks.filter_annotate_breakpoints,
        args=(
            mgd.TempInputFile("breakpoints.csv"),
            mgd.TempInputFile("library.csv"),
            [sample_id + 'N'],  # control_ids
            mgd.OutputFile(breakpoints),
            mgd.OutputFile(library),
        )
    )

    return workflow

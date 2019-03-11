function [CRseq, NCRseq]=getCRNCR(hbb,CDSindex)
    CDS = hbb.CDS;
    CRI=[];
    for i=1:2:length(CDS(CDSindex).indices)
        CRI=[CRI CDS(CDSindex).indices(i):CDS(CDSindex).indices(i+1)]; % get the coding region
        NRI=setdiff(CDS(CDSindex).indices(1):CDS(CDSindex).indices(end),CRI); % non-coding region
    end
    CRseq=hbb.Sequence(CRI);
    NCRseq=hbb.Sequence(NRI);
end
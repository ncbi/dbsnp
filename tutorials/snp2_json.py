import urllib.request
import json

class Snp2Json(object):

    api_base = 'https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/'

    def get_json(self, rs):
        
        with urllib.request.urlopen(self.api_base + str(rs)) as response:
            return(json.loads(response.read().decode()))


    def get_ss_info(self, rs_obj):
        
        if 'primary_snapshot_data' in rs_obj:
            print("\t".join(["rs", "handle", "type", "ss or RCV"]))
            self._getSsInfo(rs_obj['refsnp_id'], \
                            rs_obj['primary_snapshot_data']['support'])
            
        
    def _getSsInfo(self, rs, obj):
        for ss in obj:
            id = ss['id']
            print("\t".join([str(rs), ss['submitter_handle'],
                             id['type'], id['value']]))

        
    # def get_allele_info(self, info):
        
    #     getPlacements(info)

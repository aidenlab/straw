import json
import synapseclient
import straw
import sys

if (len(sys.argv) != 2):
   print("Usage: python {} synapse_id".format(sys.argv[0]))
   sys.exit(1)

user = "your_username"
password = "your_password"
syn = synapseclient.login(user, password)

myfile=syn.get(sys.argv[1],downloadFile=False)

batch_file_request = {
    "requestedFiles": [
        {
            "fileHandleId": myfile._file_handle['id'],
            "associateObjectId": sys.argv[1],
            "associateObjectType": "FileEntity"
        }
    ],
    "includePreSignedURLs": True,
    "includeFileHandles": False
}

result = syn.restPOST(uri='/fileHandle/batch', body=json.dumps(batch_file_request), endpoint=syn.fileHandleEndpoint)
url=result["requestedFiles"][0]["preSignedURL"]

result2 = straw.straw('NONE', url, '1', '1', 'BP', 1000000)
for i in range(len(result2[0])):
   print("{0}\t{1}\t{2}".format(result2[0][i], result2[1][i], result2[2][i]))


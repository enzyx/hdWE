from iteration import Iteration

it = Iteration(0)
print(it.getId())
print(it.getNumberOfBins())

it.generateBin(0,0,0,10)
it.generateBin(1,1,1,10)
it.generateBin(0,2,3,10)

# print(it.getNumberOfBins())
for b in it:
    print(b)
    b.generateSegment(0.1,0,0)
    b.generateSegment(0.22,1,0)
    for s in b.segments:
        print(s, s.getId(), s.getBinId())
    print(b.getProbability())

print(it.getProbability())
print(it.getNumberOfBins())
print(it.getNumberOfSegments())
print(it.getTargetNumberOfSegments())

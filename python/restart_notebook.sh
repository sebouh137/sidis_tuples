nohup jupyter-notebook --port=5678 &
sleep 5
grep token nohup.out
